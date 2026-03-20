#!/bin/bash
#SBATCH --account=small_grains
#SBATCH --partition=kgwasflow_test_elly
#SBATCH --job-name="S3"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=200GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -e

# Usage: bash 03_combine_and_filter_kmers.sh <KMER_LENGTH> <MAC> <PERCENT>
# Example: bash 03_combine_and_filter_kmers.sh 31 5 0.2

KMER_LENGTH="$1"
MAC="$2"      # e.g., 5
PERCENT="$3"  # e.g., 0.2

if [[ -z "$KMER_LENGTH" || -z "$MAC" || -z "$PERCENT" ]]; then
  echo "Usage: bash 03_combine_and_filter_kmers.sh <KMER_LENGTH> <MAC> <PERCENT>" >&2
  echo "Example: bash 03_combine_and_filter_kmers.sh 31 5 0.2" >&2
  exit 1
fi

# Generate kmers_list_paths.txt from both FT10_2/kmc and FT10_3/kmc
# Each line: full_path_to_kmers_with_strand<tab>sample_name
rm -f kmers_list_paths.txt
kmc_roots=(
  "/90daydata/small_grains/Jinglie/dataForKmerGWAS/Arabidopsis/FT10_2/kmc"
  "/90daydata/small_grains/Jinglie/dataForKmerGWAS/Arabidopsis/FT10_3/kmc"
)

for kmc_root in "${kmc_roots[@]}"; do
  [[ -d "$kmc_root" ]] || continue
  for dir in "$kmc_root"/*/; do
    [[ -d "$dir" ]] || continue
    sample="$(basename "$dir")"
    f="$kmc_root/$sample/kmers_with_strand"
    if [[ -f "$f" ]]; then
      printf "%s\t%s\n" "$f" "$sample" >> kmers_list_paths.txt
    fi
  done
done
if [[ ! -s kmers_list_paths.txt ]]; then
  echo "Error: no kmers_with_strand files found in FT10_2/kmc or FT10_3/kmc. Run step 02 first." >&2
  exit 1
fi

# Sort by sample_name (2nd column) in numeric order.
sort -t $'\t' -k2,2n kmers_list_paths.txt -o kmers_list_paths.txt

# Combine and filter
./bin/list_kmers_found_in_multiple_samples -l kmers_list_paths.txt -k "$KMER_LENGTH" --mac "$MAC" -p "$PERCENT" -o kmers_to_use

echo "Filtered k-mers list created: kmers_to_use"