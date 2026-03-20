#!/bin/bash
#SBATCH --account=small_grains
#SBATCH --job-name=fastp_kmc
#SBATCH --partition=kgwasflow_test_elly
#SBATCH --array=1-876%50              # 注意：如果样本数变了，需手动修改这里的 877
#SBATCH --cpus-per-task=16            # fastp 线程数
#SBATCH --mem=200GB
#SBATCH --time=2-00:00:00
#SBATCH --output=log/fastp_%a_%j.out
#SBATCH --error=log/fastp_%a_%j.err

# --- 1. 环境激活 ---
# source $(conda info --base)/etc/profile.d/conda.sh
conda activate kmersgwas

# --- 2. 变量设置 ---
OUT_DIR="./"
INPUT_DIR="/90daydata/small_grains/Jinglie/dataForKmerGWAS/Arabidopsis/FT10_3/fastq_out"

# 检查参数
if [ -z "$OUT_DIR" ]; then
    echo "Usage: sbatch $0 <output_directory>"
    exit 1
fi

# 创建必要的目录
mkdir -p ${OUT_DIR}/fastq ${OUT_DIR}/reports log

# --- 3. 动态生成样本列表并提取当前任务的样本名 ---
# 将你的 ls 命令结果存入 Bash 数组
SAMPLES=($(ls ${INPUT_DIR}/*.fastq | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}' | sort | uniq))

# SLURM_ARRAY_TASK_ID 是从 1 开始的，而 Bash 数组是从 0 开始的
# 所以这里需要减 1
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
SAMPLE_ID=${SAMPLES[$INDEX]}

# 打印调试信息
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing Index: $INDEX"
echo "Sample ID: $SAMPLE_ID"

# 检查样本名是否获取成功
if [ -z "$SAMPLE_ID" ]; then
    echo "Error: No sample found for index $INDEX. Is the array range correct?"
    exit 1
fi

# --- 4. 运行 fastp ---
fastp --thread $SLURM_CPUS_PER_TASK \
      --in1 ${INPUT_DIR}/${SAMPLE_ID}_1.fastq \
      --in2 ${INPUT_DIR}/${SAMPLE_ID}_2.fastq \
      --out1 ${OUT_DIR}/fastq/${SAMPLE_ID}_1.fastq \
      --out2 ${OUT_DIR}/fastq/${SAMPLE_ID}_2.fastq \
      --html ${OUT_DIR}/reports/fastp_${SAMPLE_ID}.html

echo "Sample $SAMPLE_ID finished at $(date)"