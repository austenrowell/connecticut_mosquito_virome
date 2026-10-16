#!/bin/bash
#SBATCH -J coverage_only
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -o logs/coverage_only_%j.out
#SBATCH -e logs/coverage_only_%j.err
#SBATCH -p nocona
#SBATCH --export=ALL

# Input and output files
SAMPLES_FILE="metadata.txt"
OUT_FILE="coverage_summary.tsv"

# Check for samtools
if ! command -v samtools &>/dev/null; then
  echo "✗ samtools not found in PATH. Please install samtools first." >&2
  exit 1
fi

# Write header
printf "Sample\tHost\tAvgCoverage\n" > "$OUT_FILE"

# Skip header and loop
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
      echo "[$SAMPLE] ⚠️  Missing both $SAM_GZ and $BAM, skipping." >&2
      printf "%s\t%s\tNA\n" "$SAMPLE" "$HOST" >> "$OUT_FILE"
      continue
    fi
  fi

  # Run samtools coverage → per-scaffold stats
  samtools coverage -@ 4 "$BAM" > "${SAMPLE}/coverage.tsv"

  # Compute genome-wide meanCoverage weighted by scaffold lengths
  AVG=$( awk 'NR>1 {
                   # $2 = scaffold length (LN), $5 = meanCoverage
                   sumLen += $2;
                   sumCov += $5 * $2;
               }
               END {
                 if (sumLen>0) printf "%.2f", sumCov/sumLen;
                 else print "0.00";
               }' "${SAMPLE}/coverage.tsv" )

  printf "%s\t%s\t%s\n" "$SAMPLE" "$HOST" "$AVG" >> "$OUT_FILE"
done

echo "Done. Summary written to $OUT_FILE"
