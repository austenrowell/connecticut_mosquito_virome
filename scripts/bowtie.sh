#!/bin/bash
#SBATCH -J bowtie_array
#SBATCH --array=0-20
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -o logs/bowtie_%A_%a.out
#SBATCH -e logs/bowtie_%A_%a.err
#SBATCH -p nocona
#SBATCH --export=ALL

# Get sample info
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" enriched.txt)
SAMPLE=$(echo "$LINE" | cut -f1)
HOST=$(echo "$LINE" | cut -f2)

# Define paths
HOST_INDEX="reference/${HOST}/${HOST}"
R1_FILE="${SAMPLE}/${SAMPLE}_R1_001.trimmed.fastq"
R2_FILE="${SAMPLE}/${SAMPLE}_R2_001.trimmed.fastq"
OUT_PREFIX="${SAMPLE}/${SAMPLE}_no_host"
HOST_PREFIX="${SAMPLE}/${SAMPLE}_host"
SAM_OUT="${SAMPLE}/host_reads.sam"

echo "[$(date)] Processing $SAMPLE with host index: $HOST"

bowtie2-2.4.2-sra-linux-x86_64/bowtie2 \
    -x "$HOST_INDEX" \
    -p 8 \
    -1 "$R1_FILE" \
    -2 "$R2_FILE" \
    --un-conc-gz "$OUT_PREFIX.fastq" \
    -S "$SAM_OUT"


echo "[$(date)] Done with $SAMPLE"
