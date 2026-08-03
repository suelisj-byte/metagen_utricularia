###############################################################################
## Heatmaps at ALL taxonomic ranks (class, order, family, genus)
## Figures C1 and C2.
###############################################################################
library(tidyverse); library(pheatmap); library(RColorBrewer)

data_dir  <- "/Users/suelisuarez/Desktop/Uni/SS2026/data/"
MIN_READS <- 1e6
TOP_N     <- 40

## ---- read once -------------------------------------------------------------
files <- list.files(data_dir, pattern = "^S\\d{2}\\.mpa\\.txt$", full.names = TRUE)
names(files) <- sub("\\.mpa\\.txt$", "", basename(files))
all_long <- imap_dfr(files, ~ read_tsv(.x, col_names = c("lineage","count"),
                                       show_col_types = FALSE) %>%
                       mutate(sample = .y))

## ---- metadata --------------------------------------------------------------
meta <- tribble(
  ~sample, ~pond,   ~species,
  "S01","III","Uaus",
  "S02","II","Uaus",
  "S03","I","Uaus",
  "S04","I","Uaus",
  "S05","IV","Uaus",
  "S06","VII","Uaus",
  "S07","VII","Uaus",
  "S08","VI","Uaus",
  "S09","V","Uaus",
  "S10","V","Uaus",
  "S11","II","Ubrem",
  "S12","II","Ubrem",
  "S13","II","Ubrem",
  "S14","I","Ubrem",
  "S15","I","Ubrem",
  "S16","I","Ubrem",
  "S17","I","Ubrem",
  "S18","I","Ubrem",
  "S19","I","Ubrem",
  "S20","I","Ubrem",
  "S21","IV","Ubrem",
  "S22","II","Ubrem",
  "S23","II","Ubrem",
  "S24","II","Ubrem",
  "S25","II","Ubrem",
  "S26","II","Ubrem",
  "S27","II","Ubrem",
  "S28","II","Ubrem",
  "S29","II","Ubrem",
  "S30","II","Ubrem",
  "S31","II","Ubrem",
  "S32","V","Uaus",
  "S33","V","Uaus"
) %>%
  mutate(pond=factor(pond,levels=c("I","II","III","IV","V","VI","VII")),
         species=factor(species,levels=c("Ubrem","Uaus")))

pond_cols <- c(
  "I"   = "#78bc21",  # sit 0  (chartreuse)
  "II"  = "#1976D2",  # sit 1  (blue)
  "IV"  = "#64B5F6",  # sit 1  (light blue)
  "III" = "#D16002",  # sit 2  (marmelade)
  "V"   = "#E65100",  # sit 2  (dark orange)
  "VI"  = "#FF6E00",  # sit 2  (hotorange)
  "VII" = "#FAB972"   # sit 2  (calm orange)
)

species_cols <- c(
  Ubrem = "#7A3E48",  # burgundy
  Uaus  = "#E6B800"   # flower yellow
)


ann_colours <- list(
  species = species_cols,
  pond    = pond_cols
)
## ---- function: build + draw heatmap at a given rank ------------------------
make_heatmap <- function(RANK) {
  rank_label <- c(c="Class", o="Order", f="Family", g="Genus")[RANK]
  prefix <- paste0(RANK, "__")
  
  wide <- all_long %>%
    mutate(last_elem = sub(".*\\|", "", lineage)) %>%
    filter(startsWith(last_elem, prefix)) %>%
    mutate(taxon = sub(paste0("^", prefix), "", last_elem)) %>%
    group_by(sample, taxon) %>% summarise(count=sum(count), .groups="drop") %>%
    pivot_wider(names_from=sample, values_from=count, values_fill=0) %>%
    column_to_rownames("taxon") %>% as.matrix()
  
  depth <- colSums(wide)
  wide  <- wide[, colnames(wide) %in% meta$sample[ 
    meta$sample %in% names(depth)[depth >= MIN_READS]], drop=FALSE]
  wide  <- wide[rowSums(wide) > 0, , drop=FALSE]
  
  rel <- sweep(wide, 2, colSums(wide), "/") * 100
  if (!is.na(TOP_N) && nrow(rel) > TOP_N)
    rel <- rel[names(sort(rowMeans(rel), decreasing=TRUE))[1:TOP_N], ]
  plt <- log10(rel + 0.01)
  
  ann <- meta %>% filter(sample %in% colnames(plt)) %>%
    column_to_rownames("sample") %>% select(species, pond)
  
  pheatmap(
    plt, annotation_col=ann, annotation_colors=ann_colours,
    clustering_method="average",
    color=colorRampPalette(brewer.pal(9,"YlGnBu"))(100),
    border_color=NA, fontsize_row=7,
    main=paste0(rank_label, " level (log10 % relative abundance)"),
    filename=paste0("heatmap_", rank_label, ".png"),
    width=11, height=9)
}

## ---- run for all ranks -----------------------------------------------------
for (r in c("c","o","f","g")) make_heatmap(r)
