# Analysis protocol

Bioinformatic pipeline for the bachelor's thesis "What Drives Displacement in Aquatic Carnivorous Plant? Understanding the Association Between Prey Diverity in the Species Displacement of _Utricularia bremii_ by _Utricularia australis_" (Sueli Suarez Jordan, Institute for Evolution and Biodiversity of Plants, University of Münster).

This document records every step run on the University of Münser PALMA-II HPC cluster, from raw reads to the abundance tables used for the R analyses. Some steps inclued the problems encountered and how they were solved, so that the pipeline can be reproduced or adapted.

## Softwares
- **fastp (v0.24.4)** to read quality, filter and trimm raw reads
- **Kraken2 (v2.17.1)** for host-DNA removal and taxonomic classification
- **KrakenTools ([Lu et al. 2022](https://doi.org/10.1038/s41596-022-00738-y))** to combine reports, convert to MPA, and beta diversity analysis
- **R (v4.5.1)** for statistics and figures

## Smaple naming
Cluster filenames have a descriptive prefix (e.g. ``UbremTVPfl2_S33``) with the sample code (S##). After the third step, analyses use only the sample code. For more information on the species, pond and situation, see [data/metadata.csv](https://github.com/suelisj-byte/metagen_utricularia/blob/1538483ad680ef8888b78d14d4166e54a288fff3/data/metadata.csv).

## Pipeline Overview
1. Quality filter and trim raw reads (fastp)
2. Remove _Utricularia_ host DNA against a custom Kraken2 database
3. Classify the remaining reads against the arthropod Kraken2 database ([López Clinton & van der Valk 2025](https://doi.org/10.17044/scilifelab.29666605))
4. Convert Kraken2 reports to abundance tables (KrakenTools)
5. Analyse alpha and beta diversity, and taxonomy in R (see [scripts/R/](https://github.com/suelisj-byte/metagen_utricularia/tree/289cba8a9210606f54e3eb7e39dff4c15831c5cd/scripts/R))

### 1. Retrieve raw sequencing data
The dataset (33 libraries, Illumina NovaSeq X, paired-end 2 × 150 bp) was produced for a previous master's project and provided with adapters already trimmed and demultiplexed.
```bash
cp -r /cloud/wwu1/r_agmueller/r_agmueller/Ubremii/RESEQ/ /scratch/tmp/ssuarezj/
```

### 2. Quality filtering and trimming
The first 10 bases of both reads were removed to eliminate non-random base composition at the read start. 3' ends were trimmed with a sliding window. Reads shorter than 40 bp were discarded. Per sample HTML and JSON reports were written to ``FASTP/reports/``.

For the script see [scripts/slurm/fastp.sh](https://github.com/suelisj-byte/metagen_utricularia/blob/68c7fa59f4b03f3b4843e13f5e4679447e9ff757/scripts/slurm/fastp.sh).

### 3. Host-DNA removal
#### 3.1 Build the _Utricularia_ library
Three genome assembies were combined: _U. gibba_ (NCBI GCA_002189035.1), _U. reniformis_ (NCBI GCA_009725065.1), and an unpublished draft _U. bremii_ assembly from the research group.
```bash
# U. gibba
wget "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_002189035.1/download?include_annotation_type=GENOME_FASTA"
# U. reniformis
wget "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_009725065.1/download?include_annotation_type=GENOME_FASTA"
```

The _U. bremmi_ draft assembly (``scaffolds_equalorlargerthan_1000.filtered.fasta``) is an unpublished, group-internal file and is not redistributed here. Contact the [research group](https://www.uni-muenster.de/Evolution/plantsevolbiodiv/people/index.html) for access.

**Note:** The genomic FASTA is at ``ncbi_dataset/data/<accession>/*_genomic.fna``inside it, as the NCBI datasets endpoit returns a ZIP, even when a FASTA file is requested.

```bash
mkdir kraken2_Utricularia_db
```

Taxonomy Download Option A which connects to [NCBI](https://ftp.ncbi.nlm.nih.gov)
```bash
kraken2-build --download-taxonomy \
--db kraken2_Utricularia_db
```

Taxonomy Download Option B (in case PALMA cannot reach [NCBI](https://ftp.ncbi.nlm.nih.gov))
```bash
cd kraken2_Utricularia_db && mkdir taxonomy && cd taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
cd ..
```

FASTA headers were tagged with NCBI TaxIDs (found with ``grep "<species>" names.dmp``):
| Species  | TaxID |  Header command    |
|-------|-----|-------|
| _U. bremii_ | 1896637  | sed 's/^>/>kraken:taxid|1896637|/' U_bremii.fasta > U_bremii_tax.fasta |
| _U. reniformis_   | 192314  | sed 's/^>/>kraken:taxid|192314|/' U_ren.fasta > U_ren_tax.fasta   |
| _U. gibba_   | 13748  |  sed 's/^>/>kraken:taxid|13748|/' U_gibba.fasta > U_gibba_tax.fasta  |

```bash
cat U_bremii_tax.fasta U_gibba_tax.fasta U_ren_tax.fasta > Utricularia_combined.fasta
kraken2-build --add-to-library Utricularia_combined.fasta --db kraken2_Utricularia_db
cp taxonomy/*.dmp kraken2_Utricularia_db/taxonomy/
kraken2-build --build --threads 8 --db kraken2_Utricularia_db
```

**Note:** For better results, use a _U. australis_ genomic assembly instead of _U. reniformis_ and _U. gibba_ assemblies, as host removal wil be equally efficient for both species.

#### 3.2 Run host removal
Classified plant reads are discarded, whilst unclassified reads are kept for arthropod classification.

For the script see [scripts/slurm/kraken2_host.sh](https://github.com/suelisj-byte/metagen_utricularia/blob/1538483ad680ef8888b78d14d4166e54a288fff3/scripts/slurm/kraken2_host.sh)

For per-sample values see [data/master_sequencing_summary.csv](https://github.com/suelisj-byte/metagen_utricularia/blob/1538483ad680ef8888b78d14d4166e54a288fff3/data/master_sequencing_summary.csv)

### 4. Arthropod reference database
Downloaded from [López Clinton & van der Valk (2025)](https://doi.org/10.17044/scilifelab.29666605). The database is about 1.1 TB when compressed (2,593 arthropod reference from NCBI by March 2023) and was designed for genus-level classification of Swedish arthropods.
**Note:** Such database design limits the southern-Bavarian samples in this study.

Files were fetched individually from figshare. Among all files, only ``hash.k2d.gz``(~600 GB) takes a long time to download, so a persistent session (like ``tmux`` of ``screen``) is advisable for it.

```bash
# generic pattern used for each file (hash.k2d.gz, taxo.k2d, opts.k2d,seqid2taxid.map.gz, genome_assembly_metadata.tsv):
api_url="https://api.figshare.com/v2/articles/29666605/files/<FILE_ID>"
real_url=$(curl -s "$api_url" | jq -r '.download_url')
aria2c -c -x 1 -s 1 -o <OUTPUT_NAME> "$real_url"

# Decompress hash.k2d.gz
gunzip hash.k2d.gz
```
**Note:** The decompressed hash may exhaust scratch space, so compressing the earlier host-step outputs (``clean_reads/`` and ``plant_reads/``) frees room. 

Test the database loads with ``kraken2-inspect --db /path/to/REF/Arthropod | head``

### 5. Taxonomic classification against the arthropod database
Kraken2 loads the whole hash table into RAM, so a high-memory node is required. Available partitions were checked with ``sinfo -p highmem`` and ``sinfo -o "%P %m %l %c"``, consequently ``largesmp``node (3 TB RAM) was used.

Each sample was submitted as an independent Slurm job instead of a single loop, enabling concurrent execution and freeing resources, reducing wall time from about a week to less than a day.

For the template of one job, see [scripts/slurm/kraken2_arthropod.slurm](https://github.com/suelisj-byte/metagen_utricularia/blob/1538483ad680ef8888b78d14d4166e54a288fff3/scripts/slurm/kraken2_arthropod.slurm).

Submit and monitor:
```bash
for i in sample*.slurm; do sbatch "$i"; done

# sanity check: each job writes a unique report/output
for i in sample*.slurm; do
    echo "==== $i ===="; grep "job-name\|report\|kraken_output" "$i"
done

watch -n 60 'squeue -u username'
```

### 6. Convert Kraken2 reports to abundance tables
```bash
# combined matrix with per-sample columns
combine_kreports.py --display-headers -r S*.report -o matrix_arthropod.report

# single pooled community profile (optional cross-check)
combine_kreports.py --only-combined -r S*.report -o total_arthropod.report

# per-sample MetaPhlAn-style tables (input to the R scripts)
for r in S*.report; do
    kreport2mpa.py -r "$r" -o "${r%.report}.mpa.txt"
done
```

Transfer the MPA tables to the local machine for the R analysis:
```bash
scp -J sebb00 username@palma-login:/path/to/KRAKEN2_arthropods/mpa/*.mpa.txt \
    ~/local/path/data/
```

### 7. Statistical analysis and figures
All diversity statistics and figures were produced from the ``.mpa.txt``tables in [data/mpa](https://github.com/suelisj-byte/metagen_utricularia/tree/1538483ad680ef8888b78d14d4166e54a288fff3/data/mpa). Scripts of the analyses are in [scripts/R](https://github.com/suelisj-byte/metagen_utricularia/tree/1538483ad680ef8888b78d14d4166e54a288fff3/scripts/R) and the results in . 

| Script  | Produces |
|-------|-----|
| [scripts/R/figures_crustacea.R](https://github.com/suelisj-byte/metagen_utricularia/blob/08c32311c388dbf529a5cef21d4e8f2220b273ed/scripts/R/figures_crustacea.R) | Crustacean main analyses: Shannon alpha diversity (Figure 1, Table 2 and 3); NMDS (Figure 3, NMDS stress in Table B1); PCoA (Figure B1); relative-abundance bars (Figure 4); Wilcoxon rank-sums, Kruskal-Wallis, PERMANOVA, PERMDISP, and within-Pond-I tests (Table 3) |
| [scripts/R/arthropod_ranks_analysis.R](https://github.com/suelisj-byte/metagen_utricularia/blob/d01fdaeb291c3a5ae2efd4e0af4906f183b810bb/scripts/R/arthropod_ranks_analysis.R)  | Whole-arthropod analyses: Shannon alpha diversities (all tables and figures in Appendix A); NMDS and PCoA (Figure B2) |
|  [scripts/R/crust_heatmap.R](https://github.com/suelisj-byte/metagen_utricularia/blob/d01fdaeb291c3a5ae2efd4e0af4906f183b810bb/scripts/R/crust_heatmap.R)  | Crustacean heatmaps in genus and family levels (Figure 2)  |
|  [scripts/R/arthropod_ranks_heatmap.R](https://github.com/suelisj-byte/metagen_utricularia/blob/d01fdaeb291c3a5ae2efd4e0af4906f183b810bb/scripts/R/arthropod_ranks_heatmap.R)  | Whole-arthropod heatmaps at class, order, family, and genus levels (Appendix C, Figs C1 and C2) |

Two abundance conventions exist and shouldn't be confused. Whitelisted crustacean analyses calculate proportions within crustaceans, while all-arthropod heatmaps show proportions within each rank across all arthropods. Heatmap cells are ``log10(%+0.01)``, displaying the top 40 taxa by mean relative abundance.
