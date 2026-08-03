#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=40G
#SBATCH --partition=normal
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/scratch/tmp/ssuarezj/KRAKEN2_host
#SBATCH --job-name=remove_Utricularia
#SBATCH --output=kraken_host_%j.out

module purge
source /path/to/conda.sh
conda activate kraken2
set -e

THREADS=8
DB=/path/to/kraken2_Utricularia_db
INPUT=/path/to/FASTP/trimmed
OUT=/path/to/KRAKEN2_host

mkdir -p ${OUT}/{reports,clean_reads,plant_reads,kraken_output}

for R1 in ${INPUT}/*R1.trimmed.fastq.gz; do
    sample=$(basename "$R1" _R1.trimmed.fastq.gz)
    R2=${INPUT}/${sample}_R2.trimmed.fastq.gz
    echo "Processing $sample"
    kraken2 \
        --db $DB --paired --gzip-compressed --threads $THREADS --confidence 0.1 \
        --report ${OUT}/reports/${sample}.report \
        --output ${OUT}/kraken_output/${sample}.kraken \
        --classified-out ${OUT}/plant_reads/${sample}_#.plant.fastq \
        --unclassified-out ${OUT}/clean_reads/${sample}_#.clean.fastq \
        $R1 $R2
done
