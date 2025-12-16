#!/bin/bash -l
#SBATCH --job-name=hic_download_ENA
#SBATCH --output=/scratch/grp/msc_appbio/Group2_ABCC/HiC/logs/%x_%j.out
#SBATCH --error=/scratch/grp/msc_appbio/Group2_ABCC/HiC/logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -p msc_appbio

cd /scratch/grp/msc_appbio/Group2_ABCC/HiC/raw_fastq_HiC

echo "Starting ENA FASTQ downloads..."

for srr in SRR6678863 SRR6678864 SRR6678865 SRR6678866 SRR6678867 SRR6678868 SRR6678869 SRR6678870
do
    echo "=== Getting URLs for $srr ==="

    urls=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr}&result=read_run&fields=fastq_ftp&format=tsv" | tail -n 1 | cut -f2)

    # If for some reason ENA didn't return anything, skip
    if [[ -z "$urls" ]]; then
        echo "No URLs found for $srr, skipping."
        continue
    fi

    for url in $(echo "$urls" | tr ';' '\n'); do
        echo "Downloading ftp://$url"
        wget -nc ftp://$url
    done

    echo "=== Done with $srr ==="
done

echo "All ENA downloads completed."
