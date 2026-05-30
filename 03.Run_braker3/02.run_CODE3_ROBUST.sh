#!/bin/bash
# ==============================================================================
# 代码 3 (自适应终极版): 智能 BAM 选择 + 层级降级 TSEBRA 流程
# 逻辑: 
#   1. BAM 选择优先顺: 
#      a. 自身 ID 匹配
#      b. 目录中 MapQ 高质量比对数最多的 BAM (近缘种)
#   2. 运行 BRAKER1 (RNA): 如果失败，不退出，自动标记为 RNA_FAIL。
#   3. 运行 BRAKER2 (Prot): 始终运行。
#   4. TSEBRA 合并: 
#      - 若 RNA 成功: 合并 RNA + Prot
#      - 若 RNA 失败: 仅输出 Prot 结果 (实现自动降级)
# ==============================================================================

# --- 配置 ---
SIF_IMAGE="/share/org/YZWL/yzwl_liuchao/Software/braker3.sif"
AUGUSTUS_CONFIG_BASE="./augustus_config_writable" 
GENEMARK_KEY="$HOME/.gm_key"
SAMTOOLS_EXE="$HOME/miniconda3/envs/braker3/bin/samtools"

# --- 输入 ---
TARGET_ID=${1:-"test_species_adaptive"}
GENOME_FA=${2:-"genome.fasta"}
BAM_DIR=${3:-"path_to_bam_folder"} 
PROT_SEQ=${4:-"proteins.fasta"}

THREADS=32
MAIN_WORKDIR="braker_out_${TARGET_ID}_ADAPTIVE"

# 阈值：如果最佳 BAM 的高质量 reads 少于此值，直接跳过 RNA 模式，防止 GeneMark 崩溃
MIN_BAM_READS=100000 

# ==============================================================================
# 1. 准备数据
# ==============================================================================
echo ">>> [Step 1] 准备数据..."
mkdir -p "temp_adapt_${TARGET_ID}"
mkdir -p "$MAIN_WORKDIR"
CONFIG_ABS=$(mkdir -p "$AUGUSTUS_CONFIG_BASE" && readlink -f "$AUGUSTUS_CONFIG_BASE")
KEY_ABS=$(readlink -f "$GENEMARK_KEY")

# 1.1 清洗基因组、过滤极短序列并清理非法空格
GENOME_CLEAN="temp_adapt_${TARGET_ID}/genome_cleaned.fa"
echo ">>> 正在过滤短于 3000bp 的碎片序列，并清洗序列间隙..."
# 使用 seqkit -m 3000 确保留给 GeneMark 足够长的训练序列
# 使用 -g 移除序列内部的非法空格和 gaps，防止 GeneMark 崩溃
seqkit seq -m 3000 -g "$GENOME_FA" | sed 's/[[:space:]].*$//' > "$GENOME_CLEAN"
GENOME_ABS=$(readlink -f "$GENOME_CLEAN")

# 1.2 清洗蛋白 (深度清洗 ID)
PROT_CLEAN="temp_adapt_${TARGET_ID}/proteins_cleaned.fa"
awk '/^>/ {print $1; next} {print}' "$PROT_SEQ" | sed 's/|/_/g' > "$PROT_CLEAN"
PROT_ABS=$(readlink -f "$PROT_CLEAN")

# 1.3 碎片化基因组检测与线程自适应
SCAFFOLD_COUNT=$(grep -c "^>" "$GENOME_ABS")
if [ "$SCAFFOLD_COUNT" -gt 8000 ]; then
    echo ">>> [Warning] 过滤后基因组依然高度碎片化 (包含 $SCAFFOLD_COUNT 个 scaffold)！"
    echo ">>> 自动将线程数降级为 1 (Linear Mode)。"
    THREADS=1
else
    echo ">>> 过滤后基因组连续性良好 ($SCAFFOLD_COUNT 个 scaffold)，使用多线程: 32"
    THREADS=32
fi

# ==============================================================================
# 2. 智能 BAM 选择逻辑 (Select Smart BAM)
# ==============================================================================
echo ">>> [Step 2] 正在挑选最佳 BAM..."

USE_BAM_MODE=0
SELECTED_BAM=""
BAM_ABS=""
BAI_ABS=""

if [ -d "$BAM_DIR" ]; then
    # --- 策略 A: ID 匹配 (Priority 1) ---
    MATCHED_BAM=$(find "$BAM_DIR" -name "*${TARGET_ID}*.bam" | grep -v "merged" | head -n 1)
    
    if [ ! -z "$MATCHED_BAM" ]; then
        echo ">>> [Strategy A] 找到自身 ID 匹配的 BAM: $(basename "$MATCHED_BAM")"
        SELECTED_BAM="$MATCHED_BAM"
    else
        # --- 策略 B: 最佳近缘种 (Priority 2) ---
        echo ">>> [Strategy B] 未找到自身 BAM。正在扫描目录寻找最佳近缘种 BAM..."
        echo ">>> 标准: High Quality Reads (MapQ >= 30) 数量最多"
        
        # 查找所有 BAM
        ALL_BAMS=$(find "$BAM_DIR" -name "*.bam" | grep -v "merged")
        
        MAX_HQ_READS=0
        BEST_CANDIDATE=""
        
        for bam in $ALL_BAMS; do
            # 统计 MapQ >= 30 的 reads 数，这代表高置信度比对，比单纯 reads 数更可靠
            COUNT=$("$SAMTOOLS_EXE" view -c -F 4 -q 30 "$bam")
            
            # 可选：打印进度
            # echo "  Scanning: $(basename "$bam") - HQ Reads: $COUNT"
            
            if [ "$COUNT" -gt "$MAX_HQ_READS" ]; then
                MAX_HQ_READS=$COUNT
                BEST_CANDIDATE="$bam"
            fi
        done
        
        if [ "$MAX_HQ_READS" -ge "$MIN_BAM_READS" ]; then
            echo ">>> [Strategy B] 选中最佳近缘种 BAM: $(basename "$BEST_CANDIDATE")"
            echo ">>> 高质量 Reads 数: $MAX_HQ_READS"
            SELECTED_BAM="$BEST_CANDIDATE"
        else
            echo ">>> [Warning] 即使是最佳 BAM，数据量也过低 ($MAX_HQ_READS < $MIN_BAM_READS)。"
            echo ">>> 将跳过 RNA 步骤，直接进入 Protein 模式。"
        fi
    fi
    
    # 如果选中了 BAM，建立链接和索引
    if [ ! -z "$SELECTED_BAM" ]; then
        BAM_LINK="temp_adapt_${TARGET_ID}/training.bam"
        ln -sf "$(readlink -f "$SELECTED_BAM")" "$BAM_LINK"
        
        # 检查/建立索引
        if [ -f "${SELECTED_BAM}.bai" ]; then
            ln -sf "$(readlink -f "${SELECTED_BAM}.bai")" "${BAM_LINK}.bai"
        elif [ -f "${SELECTED_BAM%.*}.bai" ]; then
            ln -sf "$(readlink -f "${SELECTED_BAM%.*}.bai")" "${BAM_LINK}.bai"
        else
            echo ">>> 正在为选中 BAM 建立索引..."
            "$SAMTOOLS_EXE" index -@ "$THREADS" "$BAM_LINK"
        fi
        
        BAM_ABS=$(readlink -f "$BAM_LINK")
        BAI_ABS=$(readlink -f "${BAM_LINK}.bai")
        USE_BAM_MODE=1
    fi
else
    echo ">>> [Warning] BAM 目录不存在。将跳过 RNA 步骤。"
fi

# ==============================================================================
# 3. 运行 BRAKER1 (RNA Mode) - 可失败
# ==============================================================================
RNA_SUCCESS=0
DIR_RNA="${MAIN_WORKDIR}/run_rna"

if [ $USE_BAM_MODE -eq 1 ]; then
    echo "----------------------------------------------------------------"
    echo ">>> [Step 3] 尝试运行 BRAKER1 (RNA模式)..."
    echo "----------------------------------------------------------------"
    DIR_RNA="${MAIN_WORKDIR}/run_rna"
    rm -rf "$DIR_RNA"   # <--- 新增这行，强制清理上次失败的残留
    mkdir -p "$DIR_RNA"

    SCRIPT_RNA="${DIR_RNA}/run_rna.sh"
    cat << EOF > "$SCRIPT_RNA"
#!/bin/bash
set -e
# Config Init
TARGET_CONFIG="/opt/Augustus/config_writable"
if [ -z "\$(ls -A \$TARGET_CONFIG)" ]; then
    cp -r /opt/Augustus/config/* \$TARGET_CONFIG/ 2>/dev/null || cp -r /usr/share/augustus/config/* \$TARGET_CONFIG/
    chmod -R u+w \$TARGET_CONFIG
fi
# 【新增】：运行前强制删除可能存在的同名旧物种配置，防止报错
rm -rf \$TARGET_CONFIG/species/${TARGET_ID}_rna

cd /output

echo "Running BRAKER1 (BAM Only)..."
# 注意: 即使是近缘种，GeneMark-ET 也比 ETP 稳定得多
/opt/BRAKER/scripts/braker.pl \\
    --genome=/genome.fa \\
    --bam=/rnaseq.bam \\
    --workingdir=/output \\
    --threads=${THREADS} \\
    --species=${TARGET_ID}_rna \\
    --gff3 \\
    --softmasking \\
    --fungus \\
    --AUGUSTUS_CONFIG_PATH=\$TARGET_CONFIG
EOF
    chmod +x "$SCRIPT_RNA"

    # 使用 set +e 允许 singularity 报错而不退出主脚本
    set +e
    singularity exec --cleanenv \
        -B "${CONFIG_ABS}:/opt/Augustus/config_writable" \
        -B "${GENOME_ABS}:/genome.fa" \
        -B "${BAM_ABS}:/rnaseq.bam" \
        -B "${BAI_ABS}:/rnaseq.bam.bai" \
        -B "${DIR_RNA}:/output" \
        -B "${SCRIPT_RNA}:/run_rna.sh" \
        -B "${KEY_ABS}:${HOME}/.gm_key" \
        "${SIF_IMAGE}" /bin/bash /run_rna.sh
    
    RET_VAL=$?
    set -e

    if [ $RET_VAL -eq 0 ]; then
        echo ">>> [Success] RNA 模式运行成功！"
        RNA_SUCCESS=1
    else
        echo ">>> [Warning] RNA 模式失败 (可能 BAM 并不适配)。"
        echo ">>> 系统将自动降级，仅依赖 Protein 结果。"
        RNA_SUCCESS=0
    fi
else
    echo ">>> [Skip] 跳过 RNA 模式 (无合适 BAM)。"
fi

# ==============================================================================
# 4. 运行 BRAKER2 (Protein Mode) - 必选
# ==============================================================================
echo "----------------------------------------------------------------"
echo ">>> [Step 4] 运行 BRAKER2 (Protein模式)..."
echo "----------------------------------------------------------------"

DIR_PROT="${MAIN_WORKDIR}/run_prot"
rm -rf "$DIR_PROT"  # <--- 新增这行，强制清理上次失败的残留
mkdir -p "$DIR_PROT"

SCRIPT_PROT="${DIR_PROT}/run_prot.sh"
cat << EOF > "$SCRIPT_PROT"
#!/bin/bash
set -e
# Config Init
TARGET_CONFIG="/opt/Augustus/config_writable"
# 【新增】：运行前强制删除可能存在的同名旧物种配置，防止报错
rm -rf \$TARGET_CONFIG/species/${TARGET_ID}_prot

cd /output

echo "Running BRAKER2 (Protein Only)..."
/opt/BRAKER/scripts/braker.pl \\
    --genome=/genome.fa \\
    --prot_seq=/proteins.fa \\
    --workingdir=/output \\
    --threads=${THREADS} \\
    --species=${TARGET_ID}_prot \\
    --gff3 \\
    --softmasking \\
    --fungus \\
    --AUGUSTUS_CONFIG_PATH=\$TARGET_CONFIG
EOF
chmod +x "$SCRIPT_PROT"

PROT_SUCCESS=0
set +e
singularity exec --cleanenv \
    -B "${CONFIG_ABS}:/opt/Augustus/config_writable" \
    -B "${GENOME_ABS}:/genome.fa" \
    -B "${PROT_ABS}:/proteins.fa" \
    -B "${DIR_PROT}:/output" \
    -B "${SCRIPT_PROT}:/run_prot.sh" \
    -B "${KEY_ABS}:${HOME}/.gm_key" \
    "${SIF_IMAGE}" /bin/bash /run_prot.sh
    
if [ $? -eq 0 ]; then
    echo ">>> [Success] Protein 模式运行成功！"
    PROT_SUCCESS=1
else
    echo ">>> [Fail] Protein 模式失败。"
    PROT_SUCCESS=0
fi
set -e

# ==============================================================================
# 5. TSEBRA 合并 (根据前两步的状态决定)
# ==============================================================================
echo "----------------------------------------------------------------"
echo ">>> [Step 5] 生成最终结果..."
echo "----------------------------------------------------------------"

FINAL_GFF="${MAIN_WORKDIR}/final_annotation.gff3"

if [ $RNA_SUCCESS -eq 1 ] && [ $PROT_SUCCESS -eq 1 ]; then
    echo ">>> 模式: 混合合并 (RNA + Protein)"
    echo ">>> 调用 TSEBRA..."
    
    SCRIPT_MERGE="${MAIN_WORKDIR}/run_merge.sh"
    cat << EOF > "$SCRIPT_MERGE"
#!/bin/bash
set -e
/opt/TSEBRA/bin/tsebra.py \\
    -g /rna/braker.gtf,/prot/braker.gtf \\
    -c /opt/TSEBRA/config/default.cfg \\
    -e /rna/hintsfile.gff,/prot/hintsfile.gff \\
    -o /output/tsebra_merged.gtf
/opt/Augustus/scripts/gtf2gff.pl < /output/tsebra_merged.gtf --out=/output/final_annotation.gff3 --gff3
EOF
    chmod +x "$SCRIPT_MERGE"
    
    singularity exec --cleanenv \
        -B "${DIR_RNA}:/rna" \
        -B "${DIR_PROT}:/prot" \
        -B "${MAIN_WORKDIR}:/output" \
        -B "${SCRIPT_MERGE}:/run_merge.sh" \
        "${SIF_IMAGE}" /bin/bash /run_merge.sh

elif [ $RNA_SUCCESS -eq 1 ]; then
    echo ">>> 模式: 仅 RNA 结果 (蛋白模式失败)"
    cp "${DIR_RNA}/braker.gff3" "$FINAL_GFF"

elif [ $PROT_SUCCESS -eq 1 ]; then
    echo ">>> 模式: 仅 Protein 结果 (RNA模式跳过或失败)"
    echo ">>> 这相当于 Augustus + Protein 模式。"
    cp "${DIR_PROT}/braker.gff3" "$FINAL_GFF"
    
else
    echo ">>> [Error] 所有模式均失败。无法生成注释。"
    exit 1
fi

# 统计
if [ -f "$FINAL_GFF" ]; then
    COUNT=$(grep -c "gene" "$FINAL_GFF")
    echo "========================================================"
    echo ">>> 流程顺利完成！"
    echo ">>> 最终结果文件: $FINAL_GFF"
    echo ">>> 预测基因总数: $COUNT"
    echo "========================================================"
else
    echo ">>> 未知错误: 结果文件未生成。"
    exit 1
fi
