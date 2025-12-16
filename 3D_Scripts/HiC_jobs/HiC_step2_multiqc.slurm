#!/bin/bash -l
#SBATCH --job-name=hic_multiqc
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/HiC/logs/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/HiC/logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH -p msc_appbio

# Load MultiQC
module load py-multiqc/1.15-gcc-13.2.0-python-3.11.6

cd /scratch/grp/msc_appbio/Group2_ABCC/HiC

# Run MultiQC on the FastQC output folder
multiqc FastQC_reports/raw -o FastQC_reports
