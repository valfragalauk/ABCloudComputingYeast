#!/bin/bash -l
#SBATCH --job-name=rnaseq_R2_download
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2/logs/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2/logs/%x_%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH -p msc_appbio

# ----------------------------------------------------------
# Working directory
# ----------------------------------------------------------
cd /scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2

# ----------------------------------------------------------
# Add SRA Toolkit to PATH
# ----------------------------------------------------------
export PATH=$PATH:/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2/SRA_tools/sratoolkit.3.2.1-ubuntu64/bin

echo "SRA Toolkit version:"
fasterq-dump --version

echo "----------------------------------------------------"
echo "Starting download"
echo "Raw FASTQ will in raw_fastq/"
echo "----------------------------------------------------"

# ----------------------------------------------------------
# SY14 R2: SRR7059705
# ----------------------------------------------------------
echo "Downloading SRR7059705"
fasterq-dump SRR7059705 \
    --split-files \
    --threads ${SLURM_CPUS_PER_TASK:-4} \
    -O raw_fastq

# ----------------------------------------------------------
# BY4742 R2: SRR7059706
# ----------------------------------------------------------
echo "Downloading SRR7059706"
fasterq-dump SRR7059706 \
    --split-files \
    --threads ${SLURM_CPUS_PER_TASK:-4} \
    -O raw_fastq

echo "----------------------------------------------------"
echo "Starting compression"
echo "----------------------------------------------------"

gzip raw_fastq/*.fastq

echo "Done"
