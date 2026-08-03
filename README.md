# metagen_utricularia

Bioinformatic pipeline and analysis code for the bachelor's thesis "What Drives Displacement in Aquatic Carnivorous Plants? Understanding the Association Between Prey Diversity in the Species Displacement of _Utricularia bremii_ by _Utricularia australis_" (Sueli Suarez Jordan, Institute for Evolution and Biodiversity of Plants, University of Münster).

The study characterises the arthropod prey captured in the traps of two aquatic bladderworts, _Utricularia bremii_ HEEK EX KÖLLIKER 1839 and _U. australis_ R. BR. 1810, using shotgun metagenomic sequencing, and it asks whether prey diversity and composition differ between the two species and among ponds. This repository contains the code and the intermediate data needed to reproduce the diversity analyses and figures from the classified read tables.

## Pipeline at a glance
```
raw reads ──▶ fastp ──▶ Kraken2 (host removal) ──▶ Kraken2 (arthropod classification) ──▶ .mpa tables ──▶ R analysis
```

Full step-by-step commands, including database construction and the HPC job scripts are in [protocol.md](https://github.com/suelisj-byte/metagen_utricularia/blob/5ee42743d5fbd624860cf22770bba524610dea0d/protocol.md).

## Repository structure
```
metagen_utricularia/
├── data/
│   ├── mpa/                          # 33 per-sample abundance tables (S01–S33.mpa.txt) — input to the R scripts
│   ├── master_sequencing_summary.csv # per-sample read counts through every pipeline stage (= thesis Table 1)
│   └── metadata.csv                  # sample → species, pond, ecological situation
├── scripts/
│   ├── R/                            # analysis and figures (R 4.5.1)
│   │   ├── figures_crustacea.R       # crustacean alpha/beta diversity, ordinations, PERMANOVA/PERMDISP, Pond I
│   │   ├── crust_heatmap.R           # crustacean genus + family heatmaps
│   │   ├── arthropod_ranks_analysis.R# whole-arthropod sensitivity analysis across ranks
│   │   └── arthropod_ranks_heatmap.R # whole-arthropod heatmaps (class/order/family/genus)
│   └── slurm/                        # HPC job scripts (PALMA-II, University of Münster)
│       ├── fastp.sh
│       ├── kraken2_host.sh
│       └── kraken2_arthropod.slurm
├── protocol.md                       # full command-level protocol, raw reads → abundance tables
└── README.md
```

## Data provided here
Because the raw sequencing reads and the reference databases are too large, they are not included in this repository. What is committed is the minimal set that reproduces the analysis:
- ``data/mpa/``: The 33 classified abundance tables that the R scripts read directly.
- ``data/master_sequencing_summary.csv``: The read counts at each stage (raw → filtered → host-removed → classified), and thus the soruce of the thesis' Table 1.
- ``data/metadata.csv``: The single source of guidence for each sample's species, pond and situation.

Not included, but how to get it:
| Not in repo | Where to get it |
| --- | --- |
| Raw reads (33 libraries) | Produced for a prior project; contact the [research group](https://www.uni-muenster.de/Evolution/plantsevolbiodiv/people/index.html) |
| _Utricularia_ host database | Built from _U. gibba_ (NCBI GCA_002189035.1), _U. reniformis_ (GCA_009725065.1), and an unpublished draft _U. bremii_ assembly (See [protocol.md](https://github.com/suelisj-byte/metagen_utricularia/blob/5ee42743d5fbd624860cf22770bba524610dea0d/protocol.md)) |
| Arhtorpod database | [López Clinton & van der Valk (2025)](https://doi.org/10.17044/scilifelab.29666605), doi: 10.17044/scilifelab.29666605 |
| Kraken2 ``.report`` and/or ``.kraken`` files | Not committed given their large sizes, but regenerable from the raw reads via [protocol.md](https://github.com/suelisj-byte/metagen_utricularia/blob/5ee42743d5fbd624860cf22770bba524610dea0d/protocol.md). |

## Reproducing the figures
As the R analysis runs from the committed ``.mpa`` tables alone, no cluster access is needed.

```r
source("scripts/R/figures_crustacea.R")       # Figures 1–4, Tables 2–3
source("scripts/R/crust_heatmap.R")            # Figure 2 heatmaps
source("scripts/R/arthropod_ranks_analysis.R") # Appendix A/B
source("scripts/R/arthropod_ranks_heatmap.R")  # Appendix C
```

## Citation
If you use this pipeline or data, please cite the thesis (Suarez Jordan, 2026) and the arthropod reference database:
> López Clinton, S. & van der Valk, T. (2025). Arthropod Kraken2 Database v1. Swedish Museum of Natural History. https://doi.org/10.17044/scilifelab.29666605.v1

## Contact
Sueli Suarez Jordan: suelisuarezjordan@gmail.com
[Magnus Wolf](https://www.uni-muenster.de/Evolution/plantsevolbiodiv/people/magnuswolf.html)
[Research Group of the Institute for Evolution and Biodiversity of Plants](https://www.uni-muenster.de/Evolution/plantsevolbiodiv/people/index.html), University of Münster.
