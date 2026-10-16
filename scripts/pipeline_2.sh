#!/bin/bash
#SBATCH -J sequential_pipeline_cov
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH -o logs/sequential_cov_%j.out
#SBATCH -e logs/sequential_cov_%j.err
#SBATCH -p nocona
#SBATCH --export=ALL


#BOWTIE2_BIN="bowtie2-2.4.2-sra-linux-x86_64/bowtie2"

# Create output dirs
#mkdir -p logs \
#         kreports \
#         kraken_outputs \
#         bracken_outputs \
#         breports \
#         krona_txt \
#         plot \
#         coverage_depth
# File to accumulate per-sample average coverage
#SUMMARY_TXT="plot/coverage_summary_2.txt"
#echo -e "Sample\tAverage_Coverage" > "$SUMMARY_TXT"

while IFS=$'\t' read -r SAMPLE HOST; do
  echo "[$(date)] --- Starting $SAMPLE (host=$HOST) ---"

#  # 1) Bowtie2 host‐removal
#  HOST_INDEX="reference/${HOST}/${HOST}"
#  R1="${SAMPLE}/${SAMPLE}_R1_001.trimmed.fastq"
#  R2="${SAMPLE}/${SAMPLE}_R2_001.trimmed.fastq"
#  OUT_PREFIX="${SAMPLE}/${SAMPLE}_no_host"
#  SAM_OUT="${SAMPLE}/host_reads.sam"

#  $BOWTIE2_BIN \
#      -x "$HOST_INDEX" \
#      -p 12 \
#      -1 "$R1" \
#      -2 "$R2" \
#      --un-conc "${OUT_PREFIX}.fastq" \
#      -S "$SAM_OUT"

#  echo "    • Bowtie2 produced: $SAM_OUT"

  # 1a) Coverage calculation
#  BAM_UNSORT="${SAMPLE}/host_reads.bam"
#  BAM_SORT="${SAMPLE}/host_reads.sorted.bam"
#  DEPTH_TXT="coverage_depth/${SAMPLE}.depth.txt"
#  PLOT_PNG="plot/${SAMPLE}_coverage.png"

  # Convert, sort, index
#  samtools view -@12 -Sb "$SAM_OUT" > "$BAM_UNSORT"
#  samtools sort -@12 "$BAM_UNSORT" -o "$BAM_SORT"
#  samtools index "$BAM_SORT"

  # Per‐base depth (includes zero‐coverage positions with -a)
#  samtools depth -a "$BAM_SORT" > "$DEPTH_TXT"

  # Compute average coverage and append to summary
#  AVG=$(awk '{sum+=$3; count++} END { if(count>0) printf("%.2f", sum/count); else print "0" }' "$DEPTH_TXT")
#  echo -e "${SAMPLE}\t${AVG}" >> "$SUMMARY_TXT"

  # 2) Kraken2 classification
  KRKN_RPT="kreports/${SAMPLE}.k2report"
  KRKN_OUT="kraken_outputs/${SAMPLE}.kraken2"
  KRAKEN_DB="/lustre/work/aurowell/kraken/viral/"
  $KRAKEN2_BIN \
    --db "$KRAKEN_DB" \
    --report "$KRKN_RPT" \
    --report-minimizer-data \
    --minimum-hit-groups 3 \
    "${OUT_PREFIX}.1.fastq" "${OUT_PREFIX}.2.fastq" \
    > "$KRKN_OUT"

  # 5) Gzip intermediates
#  gzip -v "$SAM_OUT"
  gzip -v "${OUT_PREFIX}.1.fastq"
  gzip -v "${OUT_PREFIX}.2.fastq"
#  rm "$BAM_UNSORT"

  echo "[$(date)] --- Finished $SAMPLE ---"
  echo

done < samples_1.txt

echo "All samples processed."
echo "Summary of average coverage per sample in $SUMMARY_TXT"
