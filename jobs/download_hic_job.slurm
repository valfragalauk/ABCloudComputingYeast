#!/bin/bash -l
#SBATCH --job-name=hic_download
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/HiC/raw_fastq_HiC/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/HiC/raw_fastq_HiC/%x_%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH -p msc_appbio

echo "Starting Hi-C downloads on: $(date)"

module load sra-tools/3.0.3-gcc-13.2.0

# All SRR run IDs
SRRS=(
SRR6678863  # SY6 R1
SRR6678862  # SY6 R2
SRR6678861  # BY4742 R1
SRR6678860  # BY4742 R2
SRR6678859  # SY14 R1
SRR6678858  # SY14 R2
SRR6678857  # SY13 R1
SRR6678856  # SY13 R2
)

for srr in "${SRRS[@]}"; do
    echo "============================================="
    echo "Processing $srr..."
    echo "---------------------------------------------"

    # 1. Download the .sra file
    prefetch $srr

    # 2. Convert to FASTQ (paired-end)
    fasterq-dump $srr --split-files --threads 8

    # 3. Compress the results
    gzip ${srr}_1.fastq
    gzip ${srr}_2.fastq

    echo "$srr completed on: $(date)"
done

echo "All Hi-C datasets downloaded and processed successfully!"
