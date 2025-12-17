#!/bin/bash -l
#SBATCH --job-name=rnaseq_build_ref
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2/logs/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2/logs/%x_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -p msc_appbio

# ----------------------------------------------------------
# Move to reference directory
# ----------------------------------------------------------
cd /scratch/grp/msc_appbio/Group2_ABCC/RNAseq/val_R2/reference

# ----------------------------------------------------------
# Activate conda env that has RSEM installed
# ----------------------------------------------------------
source ~/.bashrc
conda activate myenv

# ----------------------------------------------------------
# Load Bowtie2 module 
# ----------------------------------------------------------
module load bowtie2/2.5.1-gcc-13.2.0-python-3.11.6

# ----------------------------------------------------------
# Set genome + annotation filenames
# ----------------------------------------------------------
GENOME="Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
ANNOT="Saccharomyces_cerevisiae.R64-1-1.110.gtf"

echo "Using genome:     $GENOME"
echo "Using annotation: $ANNOT"

# ----------------------------------------------------------
# Build RSEM reference with Bowtie2
# ----------------------------------------------------------
echo "Starting"

rsem-prepare-reference \
  --bowtie2 \
  --gtf "$ANNOT" \
  "$GENOME" \
  yeast_rsem

echo "Reference build complete. Files with prefix: yeast_rsem"
