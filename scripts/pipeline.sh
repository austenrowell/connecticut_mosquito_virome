#!/bin/bash
#SBATCH -J sequential_pipeline
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH -o logs/sequential_%j.out
#SBATCH -e logs/sequential_%j.err
#SBATCH -p nocona
#SBATCH --export=ALL


# --- Tool & DB paths ---
BOWTIE2_BIN="bowtie2-2.4.2-sra-linux-x86_64/bowtie2"
KRAKEN2_BIN="/lustre/work/aurowell/kraken/apps/kraken2/kraken2"
BRACKEN_BIN="/lustre/work/aurowell/kraken/apps/Bracken-2.9/bracken"
KREPORT2KRONA_BIN="kreport2krona.py"
KRAKEN_DB="/lustre/work/aurowell/kraken/viral"
READ_LEN=100
LEVEL="S"
THRESHOLD=10

# --- Central output directories ---
mkdir -p logs \
         kreports \
         kraken_outputs \
         bracken_outputs \
         breports \
         krona_txt \
         plot

# --- Process each sample in sequence ---
while IFS=$'\t' read -r SAMPLE HOST; do
  echo "[$(date)] --- Starting $SAMPLE (host=$HOST) ---"

  # 1) Host‐removal with Bowtie2
  HOST_INDEX="reference/${HOST}/${HOST}"
  R1="${SAMPLE}/${SAMPLE}_R1_001.trimmed.fastq"
  R2="${SAMPLE}/${SAMPLE}_R2_001.trimmed.fastq"
  OUT_PREFIX="${SAMPLE}/${SAMPLE}_no_host"
  SAM_OUT="${SAMPLE}/host_reads.sam"

  $BOWTIE2_BIN \
      -x "$HOST_INDEX" \
      -p 12 \
      -1 "$R1" \
      -2 "$R2" \
      --un-conc "${OUT_PREFIX}.fastq" \
      -S "$SAM_OUT"

  echo "    • Bowtie2 produced: $SAM_OUT, ${OUT_PREFIX}.1.fastq, ${OUT_PREFIX}.2.fastq"

  # 2) Kraken2 classification
  KRKN_RPT="kreports/${SAMPLE}.k2report"
  KRKN_OUT="kraken_outputs/${SAMPLE}.kraken2"

  $KRAKEN2_BIN \
    --db "$KRAKEN_DB" \
    --report "$KRKN_RPT" \
    --report-minimizer-data \
    --minimum-hit-groups 3 \
    "${OUT_PREFIX}.1.fastq" "${OUT_PREFIX}.2.fastq" \
    > "$KRKN_OUT"

  if [[ ! -f "$KRKN_RPT" ]]; then
    echo "ERROR: Missing Kraken2 report for $SAMPLE → skipping to next."
    continue
  fi

  # 3) Bracken abundance estimation
  BRCK_OUT="bracken_outputs/${SAMPLE}.bracken"
  BREP="breports/${SAMPLE}.breport"

  $BRACKEN_BIN \
    -d "$KRAKEN_DB" \
    -i "$KRKN_RPT" \
    -r "$READ_LEN" \
    -l "$LEVEL" \
    -t "$THRESHOLD" \
    -o "$BRCK_OUT" \
    -w "$BREP"

  if [[ ! -f "$BREP" ]]; then
    echo "ERROR: Missing Bracken report for $SAMPLE → skipping to next."
    continue
  fi

  # 4) Krona plot generation
  KTXT="krona_txt/${SAMPLE}.b.krona.txt"
  KHTML="plot/${SAMPLE}.krona.html"

  $KREPORT2KRONA_BIN -r "$BREP" -o "$KTXT"
  ktImportText "$KTXT" -o "$KHTML"

  echo "    • Krona HTML: $KHTML"

  # 5) Gzip intermediate large files to save space
  gzip -v "$SAM_OUT"
  gzip -v "${OUT_PREFIX}.1.fastq"
  gzip -v "${OUT_PREFIX}.2.fastq"

  echo "[$(date)] --- Finished $SAMPLE ---"
  echo

done < enriched_1.txt

echo "All samples processed."

