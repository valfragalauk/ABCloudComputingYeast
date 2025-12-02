#!/bin/bash -l
#SBATCH --job-name=rnaseq_R2_fastqc_raw
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2/logs/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2/logs/%x_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH -p msc_appbio

#Working directory
cd /scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2

# Load FastQC 
module load fastqc/0.12.1-gcc-13.2.0


echo "Running FastQC on raw FASTQ files in raw_fastq/ ..."
echo "Files:"
ls raw_fastq

fastqc raw_fastq/*.fastq.gz \
       --outdir fastqc_raw \
       --threads ${SLURM_CPUS_PER_TASK:-4}

echo "FastQC on raw reads finished."
