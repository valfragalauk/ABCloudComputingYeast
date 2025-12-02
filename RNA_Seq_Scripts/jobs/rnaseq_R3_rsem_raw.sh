#!/bin/bash -l
#SBATCH --job-name=rnaseq_R3_rsem_raw
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/R_3/R_3_analysis/logs/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/R_3/R_3_analysis/logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH -p msc_appbio

# ----------------------------------------------------------
# Move into working directory
# ----------------------------------------------------------
cd /scratch/grp/msc_appbio/Group2_ABCC/RNAseq/R_3/R_3_analysis

# ----------------------------------------------------------
# Properly initialise and activate conda environment
# ----------------------------------------------------------
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate myenv

# ----------------------------------------------------------
# Load required modules
# ----------------------------------------------------------
module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6
module load samtools/1.17-gcc-13.2.0-python-3.11.6

echo "Using bowtie2 at:  $(which bowtie2)"
echo "Using samtools at: $(which samtools)"
samtools --version

# ----------------------------------------------------------
# Set the RSEM reference prefix
# ----------------------------------------------------------
REF=reference/yeast_rsem
echo "Using RSEM reference: $REF"

# ----------------------------------------------------------
# RUN RSEM ON SY14 R3  (SRR7059708)
# ----------------------------------------------------------
echo "Running RSEM on RAW SRR7059708 (SY14 R3)"

rsem-calculate-expression \
  --paired-end \
  --bowtie2 \
  --num-threads ${SLURM_CPUS_PER_TASK:-4} \
  ../SY14/fastq_raw/SRR7059708_1.fastq \
  ../SY14/fastq_raw/SRR7059708_2.fastq \
  $REF \
  rsem/SRR7059708_raw

# ----------------------------------------------------------
# RUN RSEM ON BY4742 R3  (SRR7059707)
# ----------------------------------------------------------
echo "Running RSEM on RAW SRR7059707 (BY4742 R3)"

rsem-calculate-expression \
  --paired-end \
  --bowtie2 \
  --num-threads ${SLURM_CPUS_PER_TASK:-4} \
  ../BY4272/fastq_raw/SRR7059707_1.fastq \
  ../BY4272/fastq_raw/SRR7059707_2.fastq \
  $REF \
  rsem/SRR7059707_raw

# ----------------------------------------------------------
# Move BAM alignment files (if RSEM produced them)
# ----------------------------------------------------------
echo "Moving BAM alignment files"
mv rsem/SRR7059708_raw.transcript.bam alignments/SRR7059708_raw.bam 2>/dev/null || true
mv rsem/SRR7059707_raw.transcript.bam alignments/SRR7059707_raw.bam 2>/dev/null || true

echo "RSEM alignment + quantification finished"
