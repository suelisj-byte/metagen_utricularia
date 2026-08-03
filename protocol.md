# Analysis protocol

Bioinformatic pipeline for the bachelor's thesis "What Drives Displacement in Aquatic Carnivorous Plant? Understanding the Association Between Prey Diverity in the Species Displacement of _Utricularia bremii_ by _Utricularia australis_" (Sueli Suarez Jordan, Institute for Evolution and Biodiversity of Plants, University of Münster).

This document records every step run on the University of Münser PALMA-II HPC cluster, from raw reads to the abundance tables used for the R analyses. Some steps inclued the problems encountered and how they were solved, so that the pipeline can be reproduced or adapted.

## Softwares
- **fastp (v0.24.4)** to read quality, filter and trimm raw reads
- **Kraken2 (v2.17.1)** for host-DNA removal and taxonomic classification
- **KrakenTools ([Lu et al. 2022](https://doi.org/10.1038/s41596-022-00738-y))** to combine reports, convert to MPA, and beta diversity analysis
- **R (v4.5.1)** for statistics and figures

## Pipeline Overview
1. Quality filter and trim raw reads (fastp)
2. Remove _Utricularia_ host DNA against a custom Kraken2 database
3. Classify the remaining reads against the arthropod Kraken2 database ([López Clinton & van der Valk 2025](https://doi.org/10.17044/scilifelab.29666605))
4. Convert Kraken2 reports to abundance tables (KrakenTools)
5. Analyse alpha and beta diversity, and taxonomy in R (see [scripts/R/](https://github.com/suelisj-byte/metagen_utricularia/tree/289cba8a9210606f54e3eb7e39dff4c15831c5cd/scripts/R))

