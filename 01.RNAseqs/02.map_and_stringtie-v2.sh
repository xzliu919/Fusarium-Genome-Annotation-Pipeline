#!/bin/bash
# ==============================================================================
# 跨物种 RNA-seq 映射与 StringTie 降噪流水线 (修正版)
# ==============================================================================

# --- 核心配置 ---
TARGET_ID=$1              # 参数1：目标物种名
GENOME_FA=$2              # 参数2：目标基因组绝对路径
SAMPLE_LIST_FILE=$3       # 参数3：包含样品ID的文本文件

# 路径配置
RNA_DATA_DIR="/share/org/YZWL/yzwl_liuchao/Fusarium_project/01_300Fusarium_project_analysis/00_356Fusarium_genome_files/129_Funsarium_RNA/132_fqs"

# 软件路径
STAR="/share/org/YZWL/yzwl_liuchao/miniconda3/envs/braker3/bin/STAR"
SAMTOOLS="/share/org/YZWL/yzwl_liuchao/miniconda3/envs/braker3/bin/samtools"
STRINGTIE="/share/org/YZWL/yzwl_liuchao/miniconda3/envs/braker3/bin/stringtie"

# 资源分配策略
# 注意：STAR非常吃内存。如果报错 Killed，请减小 MAX_JOBS 或 limitBAMsortRAM
MAX_JOBS=2
THREADS_PER_JOB=8

# --- 环境准备 ---
WORK_DIR="work_${TARGET_ID}"
INDEX_DIR="${WORK_DIR}/star_index"
BAM_DIR="${WORK_DIR}/bams"
mkdir -p "$INDEX_DIR" "$BAM_DIR"

# 检查输入
if [[ -z "$TARGET_ID" || -z "$GENOME_FA" || -z "$SAMPLE_LIST_FILE" ]]; then
    echo "Usage: $0 <Target_ID> <Genome_Fasta> <Sample_List>"
    exit 1
fi

if [ ! -f "$SAMPLE_LIST_FILE" ]; then
    echo "Error: 样品列表文件不存在: $SAMPLE_LIST_FILE"
    exit 1
fi

# 确保 STAR 使用的全局临时目录干净
rm -rf ./_STARtmp

echo ">>> [$(date)] 开始处理物种: ${TARGET_ID} <<<"

# 1. 构建 STAR 索引
if [ ! -f "${INDEX_DIR}/Genome" ]; then
    echo "[$(date)] 正在构建 STAR 索引..."
    # 确保索引构建也有独立的临时目录
    INDEX_TMP="${WORK_DIR}/index_STARtmp"
    rm -rf "$INDEX_TMP" 
    
    $STAR --runMode genomeGenerate --genomeDir "$INDEX_DIR" \
          --genomeFastaFiles "$GENOME_FA" \
          --genomeSAindexNbases 11 \
          --runThreadN 16 \
          --outTmpDir "$INDEX_TMP"
          
    rm -rf "$INDEX_TMP"
fi

# 2. 处理单个样本的函数
process_sample() {
    local sample="$1"
    local idx_dir="$2"
    local data_dir="$3"
    local out_dir="$4"
    local threads="$5"
    
    local final_bam="${out_dir}/${sample}.sorted.bam"
    
    # 检查是否已完成
    if [ -f "$final_bam" ]; then
        echo "[$(date)] 样本 ${sample} 已存在，跳过。"
        return 0
    fi
    
    local R1="${data_dir}/${sample}_1.clean.fq.gz"
    local R2="${data_dir}/${sample}_2.clean.fq.gz"
    
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "Error: 输入文件丢失: ${sample}" >&2
        return 1
    fi

    # 定义临时目录变量
    local star_tmp_dir="${out_dir}/${sample}_STARtmp"
    local sort_tmp_dir="${out_dir}/${sample}_sort_tmp"
    
    # 确保 sort 临时目录存在
    mkdir -p "$sort_tmp_dir"
    
    # --- 关键修改 ---
    # STAR 要求 --outTmpDir 指定的目录在运行前必须不存在
    # 所以这里必须是 rm -rf，而不是 mkdir
    rm -rf "$star_tmp_dir"
    
    echo "[$(date)] ${sample}: 运行STAR比对..."
    
    # 运行 STAR
    $STAR --genomeDir "$idx_dir" \
          --readFilesIn "$R1" "$R2" \
          --readFilesCommand zcat \
          --outFilterMismatchNmax 15 \
          --outFilterMismatchNoverLmax 0.1 \
          --outFilterMultimapNmax 20 \
          --alignIntronMax 50000 \
          --outSAMtype BAM Unsorted \
          --outTmpDir "$star_tmp_dir" \
          --runThreadN "$threads" \
          --limitBAMsortRAM 5000000000 \
          --outSAMattrRGline "ID:${sample} SM:${sample}" \
          --outFileNamePrefix "${out_dir}/${sample}_" \
          2> "${out_dir}/${sample}.star.log"
    
    local star_exit_code=$?
    
    # 无论成功失败，先清理 STAR 的临时目录（它可能残留占用空间）
    rm -rf "$star_tmp_dir"

    if [ $star_exit_code -ne 0 ]; then
        echo "Error: ${sample} STAR比对失败，请检查: ${out_dir}/${sample}.star.log" >&2
        return 1
    fi
    
    # 排序
    echo "[$(date)] ${sample}: 排序 BAM..."
    $SAMTOOLS sort -@ "$threads" -T "$sort_tmp_dir/sort" \
          -o "${out_dir}/${sample}.unsorted.bam" \
          "${out_dir}/${sample}_Aligned.out.bam"
    
    # 过滤并转为最终 BAM
    $SAMTOOLS view -@ 2 -b -q 10 \
          -o "$final_bam" \
          "${out_dir}/${sample}.unsorted.bam"
    
    if [ ! -f "$final_bam" ]; then
        echo "Error: ${sample} BAM生成失败" >&2
        return 1
    fi
    
    $SAMTOOLS index "$final_bam"
    
    # 清理中间文件
    rm -f "${out_dir}/${sample}_Aligned.out.bam" \
          "${out_dir}/${sample}.unsorted.bam" \
          "${out_dir}/${sample}_Log.out" \
          "${out_dir}/${sample}_Log.progress.out" \
          "${out_dir}/${sample}_SJ.out.tab"
    rm -rf "$sort_tmp_dir"
    
    echo "[$(date)] 完成: ${sample}"
}

export -f process_sample
export STAR SAMTOOLS

# 3. 并行执行
echo "[$(date)] 开始并行比对任务 (并发数: $MAX_JOBS)..."

# 读取样本
SAMPLES=()
while IFS= read -r line || [[ -n "$line" ]]; do
    line=$(echo "$line" | tr -d '[:space:]' | tr ',' '\n')
    for sample in $line; do
        [[ -n "$sample" ]] && SAMPLES+=("$sample")
    done
done < "$SAMPLE_LIST_FILE"

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "Error: 未找到样本 ID。"
    exit 1
fi

# 控制并发的循环
active_jobs=0
for sample in "${SAMPLES[@]}"; do
    while [ "$active_jobs" -ge "$MAX_JOBS" ]; do
        wait -n
        active_jobs=$((active_jobs - 1))
    done
    
    (
        process_sample "$sample" "$INDEX_DIR" "$RNA_DATA_DIR" "$BAM_DIR" "$THREADS_PER_JOB"
    ) &
    
    active_jobs=$((active_jobs + 1))
    echo "Started ${sample} (Active jobs: $active_jobs)"
done
wait

echo "[$(date)] 所有样本比对完成。"

# 4. 合并 BAM
echo "[$(date)] 生成合并列表..."
find "$BAM_DIR" -name "*.sorted.bam" > "${WORK_DIR}/bam_list.txt"
COUNT=$(wc -l < "${WORK_DIR}/bam_list.txt")

if [ "$COUNT" -eq 0 ]; then
    echo "Error: 没有生成有效的 BAM 文件，无法合并。"
    exit 1
fi

MERGED_BAM="${WORK_DIR}/${TARGET_ID}_merged.bam"
echo "[$(date)] 合并 $COUNT 个文件 -> $MERGED_BAM"
$SAMTOOLS merge -@ 16 -b "${WORK_DIR}/bam_list.txt" "$MERGED_BAM"
$SAMTOOLS index -@ 16 "$MERGED_BAM"

# 5. 下采样
FINAL_BAM="${WORK_DIR}/${TARGET_ID}_final_subsampled.bam"
echo "[$(date)] 下采样 (20%)..."
$SAMTOOLS view -@ 16 -b -s 0.2 "$MERGED_BAM" > "$FINAL_BAM"
$SAMTOOLS index -@ 16 "$FINAL_BAM"

# 6. StringTie
OUT_GTF="${TARGET_ID}_evidence_rnaseq.gtf"
echo "[$(date)] StringTie 组装 -> $OUT_GTF"
$STRINGTIE "$FINAL_BAM" -o "$OUT_GTF" -p 16 -m 200 -c 3.0

# 清理大文件 (保留最终用于 StringTie 的文件)
# rm -f "$MERGED_BAM" "${MERGED_BAM}.bai"

echo ">>> [$(date)] 流程结束 <<<"