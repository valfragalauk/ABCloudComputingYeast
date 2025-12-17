#!/bin/bash -l
#SBATCH --job-name=hic_step2_qc_trim
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/HiC/logs/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/HiC/logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH -p msc_appbio

module load fastqc

echo "FastQC on raw Hi-C reads "

#directory
cd /scratch/grp/msc_appbio/Group2_ABCC/HiC

mkdir -p FastQC_reports/raw

#Fastqc
fastqc raw_fastq_HiC/*.fastq.gz -o FastQC_reports/raw

echo "FastQC complete "
