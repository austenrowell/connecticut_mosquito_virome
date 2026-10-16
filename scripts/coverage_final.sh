#!/bin/bash
#SBATCH -J coverage_only
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -o logs/coverage_only_%j.out
#SBATCH -e logs/coverage_only_%j.err
#SBATCH -p nocona
#SBATCH --export=ALL

SAMPLES_FILE="metadata.txt"
OUT_FILE="coverage_summary.tsv"

mkdir -p logs

if ! command -v samtools &>/dev/null; then
  echo "✗ samtools not found in PATH. Please install/load samtools first." >&2
  exit 1
fi

printf "Sample\tHost\tAvgCoverage\n" > "$OUT_FILE"

tail -n +2 "$SAMPLES_FILE" | while IFS=$'\t' read -r SAMPLE HOST; do
  BAM="${SAMPLE}/host_reads.sorted.bam"
  SAM_GZ="${SAMPLE}/host_reads.sam.gz"

  if [ ! -f "$BAM" ]; then
    if [ -f "$SAM_GZ" ]; then
      echo "[$SAMPLE] Converting SAM → sorted BAM…" >&2

      samtools view -@ 4 -bS "$SAM_GZ" \
        | samtools sort -@ 4 -o "$BAM"

      samtools index "$BAM"
    else
      echo "[$SAMPLE] ⚠️ Missing both $SAM_GZ and $BAM, skipping." >&2
      printf "%s\t%s\tNA\n" "$SAMPLE" "$HOST" >> "$OUT_FILE"
      continue
    fi
  fi

  samtools coverage -@ 4 "$BAM" > "${SAMPLE}/coverage.tsv"

  AVG=$(awk 'NR>1 {
                 # samtools coverage columns:
                 # $2 = startpos
                 # $3 = endpos
                 # $7 = meandepth
                 len = $3 - $2 + 1
                 sumLen += len
                 sumDepth += $7 * len
             }
             END {
                 if (sumLen > 0) printf "%.6f", sumDepth / sumLen
                 else print "0.000000"
             }' "${SAMPLE}/coverage.tsv")

  printf "%s\t%s\t%s\n" "$SAMPLE" "$HOST" "$AVG" >> "$OUT_FILE"
done

