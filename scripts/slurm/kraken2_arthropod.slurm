#!/bin/bash
#SBATCH --partition=largesmp ##zen3 works but takes a lot more time
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=kraken_S00
#SBATCH --output=/path/to/KRAKEN2_arthropods/jobid/S00_%j.out

set -e
module purge
source /path/to/conda.sh
conda activate kraken2_env

sample=S00

echo "Starting $(date)"
echo "Sample: ${sample}"
hostname

kraken2 \
--db /path/to/REF/Arthropod \
--memory-mapping \
--paired \
--threads 32 \
--gzip-compressed \
--report /path/to/KRAKEN2_arthropods/reports/${sample}.report \
--output /path/to/KRAKEN2_arthropods/kraken_output/${sample}.kraken \
/path/to/KRAKEN2_host/clean_reads/*_${sample}__1.clean.fastq.gz \
/path/to/tmp/ssuarezj/KRAKEN2_host/clean_reads/*_${sample}__2.clean.fastq.gz

echo "Finished $(date)"
