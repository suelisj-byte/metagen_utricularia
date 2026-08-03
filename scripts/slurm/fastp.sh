#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=20G
#SBATCH --partition=normal
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/scratch/tmp/ssuarezj/FASTP/
#SBATCH --job-name=fastp_utricularia
#SBATCH --output=fastp_%j.out
#SBATCH --mail-user=YOURMAIL@uni-muenster.de
#SBATCH --mail-type=ALL

module purge
source /path/to/conda.sh
conda activate FASTPenv
set -e

mkdir -p trimmed reports
cd /path/to/RAW/

for file in *_R1_001.fastq.gz; do
    sample=$(basename "$file" _R1_001.fastq.gz)
    r1="${sample}_R1_001.fastq.gz"
    r2="${sample}_R2_001.fastq.gz"
    fastp \
        -i "$r1" -I "$r2" \
        -o "trimmed/${sample}_R1.trimmed.fastq.gz" \
        -O "trimmed/${sample}_R2.trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --trim_front1 10 --trim_front2 10 \
        --cut_tail --cut_mean_quality 25 \
        --qualified_quality_phred 20 \
        --length_required 40 \
        --thread 8 \
        --html "reports/${sample}.fastp.html" \
        --json "reports/${sample}.fastp.json"
done
