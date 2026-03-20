#!/bin/bash
#SBATCH --account=small_grains
#SBATCH --job-name=kmc_array
#SBATCH --partition=kgwasflow_test_elly
#SBATCH --array=1-877%50              # !!! 修改：1-样本总数，%20表示并行数
#SBATCH --cpus-per-task=48            # 对应原来的 THREADS
#SBATCH --mem=200GB                   # 每个任务的内存
#SBATCH --time=2-00:00:00
#SBATCH --output=log/kmc_%a_%j.out
#SBATCH --error=log/kmc_%a_%j.err

# --- 配置参数 ---
KMER_LENGTH=31             
SAMPLES_TSV="/90daydata/small_grains/Jinglie/dataForKmerGWAS/Arabidopsis/FT10_3/FT10_3_samples.tsv"
KMC_BIN="/project/small_grains/Jinlie.Zhou/dataForKmerGWAS/Arabidopsis/FT10/run_gwas_el_2/external_programs/kmc_v3"
STRAND_BIN="/project/small_grains/Jinlie.Zhou/dataForKmerGWAS/Arabidopsis/FT10/run_gwas_el_2/bin/kmers_add_strand_information"

# --- 脚本逻辑 ---
set -e
date

# 1. 确保日志目录存在
mkdir -p log

# 2. 获取当前阵列索引对应的唯一样本名
# 这里通过 awk 提取第一列(sample_name)，去重，然后按行号选取
SAMPLE_NAME=$(awk -F'\t' 'NR>1 {print $1}' "$SAMPLES_TSV" | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "$SAMPLE_NAME" ]]; then
    echo "Error: 第 ${SLURM_ARRAY_TASK_ID} 行没有找到样本名。"
    exit 1
fi

echo "正在处理任务 ID ${SLURM_ARRAY_TASK_ID}, 样本名: $SAMPLE_NAME"

# 3. 创建输出目录
OUTDIR="kmc/$SAMPLE_NAME"
mkdir -p "$OUTDIR"

# 4. 从 TSV 中提取该样本对应的所有 Fastq 文件路径
# 逻辑：匹配第一列，打印第3列(fq1)和第4列(fq2)，自动去除空格并忽略空值
awk -F'\t' -v sample="$SAMPLE_NAME" '
  NR>1 && $1==sample {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $3);
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $4);
    if ($3 != "") print $3;
    if ($4 != "") print $4
  }
' "$SAMPLES_TSV" > "$OUTDIR/input_files.txt"

if [[ ! -s "$OUTDIR/input_files.txt" ]]; then
  echo "Error: 在 $SAMPLES_TSV 中未找到样本 '$SAMPLE_NAME' 的路径" >&2
  exit 1
fi

# 5. 运行 KMC (Canonized)
echo "Running KMC Canonized..."
"$KMC_BIN" -t"$SLURM_CPUS_PER_TASK" -k"$KMER_LENGTH" -ci2 \
    @"$OUTDIR/input_files.txt" \
    "$OUTDIR/output_kmc_canon" \
    "$OUTDIR" \
    1> "$OUTDIR/kmc_canon.1" 2> "$OUTDIR/kmc_canon.2"

# 6. 运行 KMC (Non-canonized)
echo "Running KMC All..."
"$KMC_BIN" -t"$SLURM_CPUS_PER_TASK" -k"$KMER_LENGTH" -ci0 -b \
    @"$OUTDIR/input_files.txt" \
    "$OUTDIR/output_kmc_all" \
    "$OUTDIR" \
    1> "$OUTDIR/kmc_all.1" 2> "$OUTDIR/kmc_all.2"

# 7. 合并链信息 (Combine strand information)
echo "Adding strand information..."
"$STRAND_BIN" \
    -c "$OUTDIR/output_kmc_canon" \
    -n "$OUTDIR/output_kmc_all" \
    -k "$KMER_LENGTH" \
    -o "$OUTDIR/kmers_with_strand"

# 8. 清理中间文件 (可选，如果需要节省空间可以取消注释)
# rm -f "$OUTDIR"/*.kmc*

echo "Sample $SAMPLE_NAME completed at $(date)"