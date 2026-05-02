#!/bin/bash

# 定义输入文件
FASTA=$1
GTF=$2
PREFIX=$3

echo "开始处理..."

# 使用 gffread 一步到位完成所有操作：
# -E : 暴露/显示警告信息
# -T : (如果输入是gff输出gtf，不加-T默认输出GFF3)
# -g : 指定参考基因组 fasta
# -o : 输出转换后的标准 GFF3 文件
# -x : 输出 CDS 核酸序列
# -y : 输出 Protein 氨基酸序列
gffread "$GTF" -g "$FASTA" -o "${PREFIX}_standard.gff3" -x "${PREFIX}_cds.fa" -y "${PREFIX}_pep.fa"

echo "处理完成！生成了以下文件："
echo "1. 标准GFF3文件: ${PREFIX}_standard.gff3"
echo "2. CDS序列文件: ${PREFIX}_cds.fa"
echo "3. 蛋白序列文件: ${PREFIX}_pep.fa"