#!/bin/bash -l
#SBATCH --job-name=rnaseq_R2_rsem_raw
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/R_2/logs/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/R_2/logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH -p msc_appbio

# ----------------------------------------------------------
# Move into working directory
# ----------------------------------------------------------
cd /scratch/grp/msc_appbio/Group2_ABCC/RNAseq/R_2

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
# RUN RSEM ON SY14 R2  (SRR7059705)
# ----------------------------------------------------------
echo "Running RSEM on RAW SRR7059705 (SY14 R2)"

rsem-calculate-expression \
  --paired-end \
  --bowtie2 \
  --num-threads ${SLURM_CPUS_PER_TASK:-4} \
  raw_fastq/SRR7059705_1.fastq \
  raw_fastq/SRR7059705_2.fastq \
  $REF \
  rsem/SRR7059705_raw

# ----------------------------------------------------------
# RUN RSEM ON BY4742 R2  (SRR7059706)
# ----------------------------------------------------------
echo "Running RSEM on RAW SRR7059706 (BY4742 R2)"

rsem-calculate-expression \
  --paired-end \
  --bowtie2 \
  --num-threads ${SLURM_CPUS_PER_TASK:-4} \
  raw_fastq/SRR7059706_1.fastq \
  raw_fastq/SRR7059706_2.fastq \
  $REF \
  rsem/SRR7059706_raw

# ----------------------------------------------------------
# Move BAM alignment files (if RSEM produced them)
# ----------------------------------------------------------
echo "BAM alignment files"
mv rsem/SRR7059705_raw.transcript.bam alignments/SRR7059705_raw.bam 2>/dev/null || true
mv rsem/SRR7059706_raw.transcript.bam alignments/SRR7059706_raw.bam 2>/dev/null || true

echo "RSEM alignment + quantification finished"
