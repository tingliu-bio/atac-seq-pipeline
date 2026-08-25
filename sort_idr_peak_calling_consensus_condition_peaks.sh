#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 2
#SBATCH --mem=61440M
#SBATCH -t 0-10:00 # Runtime in D-HH:MM
#SBATCH -J sort_idr

work_path=$SLURM_SUBMIT_DIR
cd ${work_path}
#set -euo pipefail

peak_list="ATL43_filtered_narrowPeak_list.txt"
outdir="IDR_results"
#mkdir -p "$outdir"

module load samtools/1.20
module load bowtie2/2.4.5
module load MACS/2.2.7.1
module load sambamba/0.7.0
module load bedtools/2.27.1

if [ ! -d  $work_path/IDR_results ]
then mkdir -p $work_path/IDR_results
fi

#conda activate DeepL_python

# 1) Sort all peak files
while read -r peak_file; do
    [[ -z "$peak_file" ]] && continue

    dir=$(dirname "$peak_file")
    base=$(basename "$peak_file")

    # remove suffix
    prefix=${base%_peaks.filtered.narrowPeak}

    sorted_peak="${outdir}/${prefix}.sorted.narrowPeak"

    sort -k8,8nr "$peak_file" > "$sorted_peak"

    echo "$prefix sort done"
done < "$peak_list"


# 2) Get condition names by removing replicate suffix _1 / _2
cut -f1 "$peak_list" | while read -r peak_file; do
    base=$(basename "$peak_file")
    prefix=${base%_peaks.filtered.narrowPeak}
    condition=$(echo "$prefix" | sed -E 's/_[12]$//')
    echo "$condition"
done | sort -u > "${outdir}/conditions.txt"


# 3) Run IDR for each condition
while read -r condition; do
    rep1="${outdir}/${condition}_1.sorted.narrowPeak"
    rep2="${outdir}/${condition}_2.sorted.narrowPeak"

    if [[ ! -f "$rep1" || ! -f "$rep2" ]]; then
        echo "WARNING: missing replicate for $condition, skip"
        continue
    fi

    echo "Running IDR for $condition"

    idr \
      --samples "$rep1" "$rep2" \
      --input-file-type narrowPeak \
      --rank p.value \
      --output-file "${outdir}/idr_${condition}_output.txt" \
      --output-file-type narrowPeak \
      --soft-idr-threshold 0.05 \
      --plot

    echo "$condition IDR done"

done < "${outdir}/conditions.txt"


# 4) Optional but recommended: explicitly keep IDR <= 0.05
# column 12 = -log10(IDR), so IDR <= 0.05 means column12 >= 1.30103
for f in ${outdir}/idr_*_output.txt; do
    awk '$12 >= 1.30103' "$f" > "${f%.txt}.filtered.narrowPeak"
done


# 5) Merge all condition-level IDR peaks
cat ${outdir}/idr_*_output.filtered.narrowPeak \
  | cut -f1-3 \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "${outdir}/consensus_IDR_peaks.bed"

echo "All done: ${outdir}/consensus_IDR_peaks.bed"
