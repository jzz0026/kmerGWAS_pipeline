# 1. 加载 Bowtie 模块
module load bowtie/1.3.1

# 2. 定义参考基因组索引路径
INDEX="/project/small_grains/Jinlie.Zhou/refGenome/Arabidopsis_thaliana/TAIR10/TAIR10_index"

# 3. 从原始结果中提取 P < 0.001 的 k-mer 及其 P 值
# 这一步会自动处理你的 rs 列 (序列_编号)，只保留序列部分
echo "Extracting significant k-mers..."
zcat output/phenotype_value.assoc.txt.gz | awk 'NR>1 && $9 < 0.001 {
    split($2, a, "_"); 
    print a[1], $9 
}' > kmers_to_map.txt

# 4. 转换为 FASTA 格式
awk '{print ">kmer_"NR"\n"$1}' kmers_to_map.txt > kmers_to_map.fa

# 5. 运行 Bowtie 比对 (遵循论文逻辑方案)
# -a 报告所有位置, --best --strata 保证位置是最优的一层
echo "Mapping k-mers to TAIR10..."
bowtie -f -v 3 -a --best --strata $INDEX kmers_to_map.fa -S kmers_mapped.sam

# 提取所有非头部的行, 统计每个 k-mer ID (第一列 $1) 出现的次数, 只打印那些只出现过 1 次的行
grep -v "@" kmers_mapped.sam | awk '{count[$1]++; line[$1]=$0} END {for (id in count) if (count[id] == 1) print line[id]}' > kmers_mapped_unique.sam

# 6. 解析比对结果坐标
# 提取列：k-mer编号, 染色体, 坐标位置
grep -v "@" kmers_mapped_unique.sam | awk '{print $1, $3, $4}' > kmer_coords.txt

# 7. 合并坐标与 P 值，生成绘图输入文件
# 结果包含五列：Chromosome, Position, P-value, kmer_id, sequence
echo "Merging coordinates with P-values and adding k-mer ID + sequence..."
awk 'NR==FNR {seq["kmer_"NR]=$1; p["kmer_"NR]=$2; next} {if(p[$1]) print $2, $3, p[$1], $1, seq[$1]}' \
    kmers_to_map.txt kmer_coords.txt > manhattan_input.txt

echo "Done! File 'manhattan_input.txt' is ready for R."