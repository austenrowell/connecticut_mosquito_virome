#!/bin/bash
#SBATCH -J trimal
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH -p nocona
#SBATCH -o logs/trimal.%A_%a.out
#SBATCH -e logs/trimal.%A_%a.err
#SBATCH --export=ALL
#SBATCH --mail-user=austen.rowell@ttu.edu
#SBATCH --mail-type=ALL


# Path to adapters
ADAPTERS=../adapters.fasta

# Updated sample list
SAMPLES=(
    Sample_10906
    Sample_11261
    Sample_11273
    Sample_12147
    Sample_2242
    Sample_4323
    Sample_440
    Sample_6271
    Sample_66
    Sample_6711
    Sample_7800
    Sample_7836
    Sample_8292
    Sample_9171
    Sample_9208
    Sample_E10906
    Sample_E11273
    Sample_E12147
    Sample_E12717
    Sample_E13649
    Sample_E12664
    Sample_E1488
    Sample_E2242
    Sample_E2388
    Sample_E3888
    Sample_E4265
    Sample_E440
    Sample_E5422
    Sample_E66
    Sample_E6711
    Sample_E6915
    Sample_E7231
    Sample_E7800
    Sample_E8244
    Sample_E9208
    Sample_E9621
    Sample_E9906
    Sample_EAal-10
    Sample_EAal-11
    Sample_EAal-12
    Sample_EAal-13
    Sample_EAal-2
    Sample_EAal-3
    Sample_EAal-4
    Sample_EAal-5
    Sample_EAal-6
    Sample_EAal-7
    Sample_EAal-8
    Sample_EAal-9
    Sample_16222
    Sample_E16222
    Sample_E8933
)

# Run Trimmomatic for each sample
echo "Starting batch trimming..."
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing $SAMPLE..."
    cd "$SAMPLE" || { echo "Failed to enter directory $SAMPLE"; continue; }

    # Dynamically get the R1 and R2 files (assuming only one pair per folder)
    R1_FILE=$(ls *R1_001.fastq.gz 2>/dev/null | head -n1)
    R2_FILE=$(ls *R2_001.fastq.gz 2>/dev/null | head -n1)

    if [[ -z "$R1_FILE" || -z "$R2_FILE" ]]; then
        echo "Missing R1 or R2 file in $SAMPLE"
        cd ..
        continue
    fi

    # Output file name base
    BASE_NAME="${R1_FILE%%_R1_001.fastq.gz}"

    trimmomatic PE -threads $SLURM_CPUS_PER_TASK \
        "$R1_FILE" "$R2_FILE" \
        "${BASE_NAME}_R1_001.trimmed.fastq" "${BASE_NAME}_R1_001.untrimmed.fastq" \
        "${BASE_NAME}_R2_001.trimmed.fastq" "${BASE_NAME}_R2_001.untrimmed.fastq" \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 HEADCROP:8

    cd ..
done

echo "All trimming jobs completed."

