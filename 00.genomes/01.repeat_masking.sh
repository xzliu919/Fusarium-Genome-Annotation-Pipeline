#!/bin/bash
# 设定变量
# 1. 加载 Conda 的基础配置 (根据你的路径)
source /share/org/YZWL/yzwl_liuchao/miniconda3/bin/activate

# 2. 激活 braker3 环境
conda activate braker3

GENOME=$1
SPECIES_NAME=$2
CPU=16

# 0. 设置工作目录和路径
WORKDIR=$(pwd)
echo "当前工作目录: $WORKDIR"
echo "基因组文件: $GENOME"
echo "物种名称: $SPECIES_NAME"

# 1.1 构建重复序列数据库
echo "步骤 1.1: 构建重复序列数据库..."
/share/org/YZWL/yzwl_liuchao/miniconda3/envs/braker3/bin/BuildDatabase -name ${SPECIES_NAME}_db $GENOME

# 检查数据库是否创建成功
if [ ! -f "${SPECIES_NAME}_db.nhr" ]; then
    echo "错误: 数据库创建失败"
    exit 1
fi

# 1.2 预测重复序列家族
echo "步骤 1.2: 预测重复序列家族 (可能需要较长时间)..."
/share/org/YZWL/yzwl_liuchao/miniconda3/envs/braker3/bin/RepeatModeler -database ${SPECIES_NAME}_db -threads $CPU -LTRStruct

# 检查输出文件
if [ ! -f "${SPECIES_NAME}_db-families.fa" ]; then
    echo "错误: 重复序列家族预测失败"
    exit 1
fi

echo "重复序列家族预测完成，发现 $(grep -c '^>' ${SPECIES_NAME}_db-families.fa) 个家族"

# 1.3 对基因组进行软屏蔽
echo "步骤 1.3: 对基因组进行软屏蔽..."

# 先检查 RepeatMasker 的可用选项
echo "检查 RepeatMasker 版本和选项..."
/share/org/YZWL/yzwl_liuchao/miniconda3/envs/braker3/bin/RepeatMasker -h 2>&1 | grep -A5 "OPTIONS" | head -10

# 创建输出目录
MASKING_OUTDIR="repeat_masking_out_${SPECIES_NAME}"
mkdir -p $MASKING_OUTDIR

# 对于真菌基因组，建议使用 -species fungi 或使用自定义库
# 方法1: 使用自定义库（从 RepeatModeler 生成）
echo "使用自定义重复序列库进行屏蔽..."
/share/org/YZWL/yzwl_liuchao/miniconda3/envs/braker3/bin/RepeatMasker \
    -lib ${SPECIES_NAME}_db-families.fa \
    -pa $(($CPU/4))  -xsmall \
    -gff \
    -dir $MASKING_OUTDIR \
    $GENOME

# 如果上述方法失败，尝试方法2: 使用物种特异性库
if [ ! -f "$MASKING_OUTDIR/${GENOME}.masked" ]; then
    echo "自定义库方法可能失败，尝试使用真菌库..."
    /share/org/YZWL/yzwl_liuchao/miniconda3/envs/braker3/bin/RepeatMasker \
        -species fungi \
        -pa $(($CPU/4)) \
        -xsmall \
        -gff \
        -dir $MASKING_OUTDIR \
        $GENOME
fi

# 输出文件
MASKED_GENOME="$MASKING_OUTDIR/${GENOME}.masked"

# 检查结果
if [ -f "$MASKED_GENOME" ]; then
    echo "重复序列屏蔽完成!"
    echo "原始基因组大小: $(grep -v '^>' $GENOME | tr -d '\n' | wc -c) bp"
    echo "屏蔽后基因组: $MASKED_GENOME"
    
    # 统计屏蔽率
    echo "计算重复序列比例..."
    total_bases=$(grep -v '^>' $MASKED_GENOME | tr -d '\n' | wc -c)
    lowercase_bases=$(grep -v '^>' $MASKED_GENOME | tr -d '\n' | grep -o '[a-z]' | wc -l)
    
    if [ $total_bases -gt 0 ]; then
        repeat_percent=$(echo "scale=2; $lowercase_bases * 100 / $total_bases" | bc)
        echo "重复序列比例: $repeat_percent%"
    fi
else
    echo "警告: 未找到屏蔽后的基因组文件"
    ls -la $MASKING_OUTDIR/
fi

echo "处理完成!"
