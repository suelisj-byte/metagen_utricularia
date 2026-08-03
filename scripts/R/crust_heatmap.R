###############################################################################
## Freshwater crustacean prey heatmaps: family and genus level
## Figure 2
###############################################################################
packages <- c("tidyverse", "pheatmap", "RColorBrewer", "grid")
to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
lapply(packages, library, character.only = TRUE)

data_dir  <- "~/path/to/mpa_data/"
MIN_READS <- 1e6

## ---- FRESHWATER WHITELIST -------------
whitelist <- tribble(
  ~family,            ~group,
  "Daphniidae","Cladocera","Chydoridae","Cladocera","Bosminidae","Cladocera",
  "Sididae","Cladocera","Macrothricidae","Cladocera",
  "Cyclopidae","Copepoda","Cyclopettidae","Copepoda","Diaptomidae","Copepoda",
  "Temoridae","Copepoda","Canthocamptidae","Copepoda",
  "Cyprididae","Ostracoda","Candonidae","Ostracoda","Darwinulidae","Ostracoda",
  "Gammaridae","Amphipoda","Crangonyctidae","Amphipoda","Hyalellidae","Amphipoda"
)

## ---- READ .mpa FILES --------------------------------------------------
files <- list.files(data_dir, pattern = "^S\\d{2}\\.mpa\\.txt$", full.names = TRUE)
names(files) <- sub("\\.mpa\\.txt$", "", basename(files))
all_long <- imap_dfr(files, ~ read_tsv(.x, col_names = c("lineage","count"),
                                       show_col_types = FALSE) %>%
                       mutate(sample = .y))

## ---- METADATA + COLOURS ----------------------------------------------------
meta <- tribble(
  ~sample, ~pond,   ~species,
  "S01","III","Uaus","S02","II","Uaus","S03","I","Uaus","S04","I","Uaus",
  "S05","IV","Uaus","S06","VII","Uaus","S07","VII","Uaus","S08","VI","Uaus",
  "S09","V","Uaus","S10","V","Uaus","S11","II","Ubrem","S12","II","Ubrem",
  "S13","II","Ubrem","S14","I","Ubrem","S15","I","Ubrem","S16","I","Ubrem",
  "S17","I","Ubrem","S18","I","Ubrem","S19","I","Ubrem","S20","I","Ubrem",
  "S21","IV","Ubrem","S22","II","Ubrem","S23","II","Ubrem","S24","II","Ubrem",
  "S25","II","Ubrem","S26","II","Ubrem","S27","II","Ubrem","S28","II","Ubrem",
  "S29","II","Ubrem","S30","II","Ubrem","S31","II","Ubrem","S32","V","Uaus",
  "S33","V","Uaus"
) %>%
  mutate(
    pond=factor(pond,levels=c("I","II","III","IV","V","VI","VII")),
    species=factor(species,levels=c("Ubrem","Uaus")),
    situation = case_when(
      pond == "I"                 ~ "0_cooccur",
      pond == "II"                ~ "1_bremii",
      pond %in% c("V","VI","VII") ~ "2_australis",
      TRUE                        ~ "other"          # III, IV -> descriptive only
    ) %>% factor(levels = c("0_cooccur","1_bremii","2_australis","other"))
  )

pond_cols <- c(
  "I"   = "#78bc21",  # sit 0  (chartreuse)
  "II"  = "#1976D2",  # sit 1  (blue)
  "IV"  = "#64B5F6",  # sit 1  (light blue)
  "III" = "#D16002",  # sit 2  (marmelade)
  "V"   = "#E65100",  # sit 2  (dark orange)
  "VI"  = "#FF6E00",  # sit 2  (hot orange)
  "VII" = "#FAB972"   # sit 2  (calm orange)
)

species_cols <- c(Ubrem="#7A3E48", Uaus="#E6B800")
group_cols   <- c(Cladocera="#3C3D3C", Copepoda="#A2A5A3",
                  Ostracoda="#CFD0CF", Amphipoda="#FCFCFC")

## ===========================================================================
## FUNCTION: build, draw, and save a freshwater-crustacean heatmap at a rank
## ===========================================================================
make_crust_heatmap <- function(LEVEL = "f") {
  
  rank_label <- c(f = "family", g = "genus")[LEVEL]
  prefix     <- paste0(LEVEL, "__")
  
  ## 1. rows at chosen rank, carry family for whitelist join -----------------
  tax_long <- all_long %>%
    mutate(last_elem = sub(".*\\|", "", lineage)) %>%
    filter(startsWith(last_elem, prefix)) %>%
    mutate(
      taxon  = sub(paste0("^", prefix), "", last_elem),
      family = ifelse(grepl("\\|f__", lineage),
                      sub(".*\\|f__([^|]+).*", "\\1", lineage), NA)
    ) %>%
    select(sample, family, taxon, count)
  
  ## 2. keep only whitelisted freshwater families ---------------------------
  fw <- tax_long %>% inner_join(whitelist, by = "family")
  if (nrow(fw) == 0) stop("No whitelisted taxa at level ", LEVEL)
  
  ## 3. wide matrix (taxon x sample) ----------------------------------------
  fw_wide <- fw %>%
    group_by(sample, taxon) %>% summarise(count = sum(count), .groups="drop") %>%
    pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
    column_to_rownames("taxon") %>% as.matrix()
  
  ## 4. drop failed libraries (total reads) ---------------------------------
  tot_depth <- all_long %>% group_by(sample) %>% summarise(d = sum(count))
  keep <- tot_depth$sample[tot_depth$d >= MIN_READS]
  fw_wide <- fw_wide[, intersect(colnames(fw_wide), keep), drop = FALSE]
  fw_wide <- fw_wide[rowSums(fw_wide) > 0, , drop = FALSE]
  cat("[", rank_label, "] taxa:", nrow(fw_wide),
      "| samples:", ncol(fw_wide), "\n")
  
  ## 5. relative abundance within crustacea, log10 --------------------------
  fw_rel <- sweep(fw_wide, 2, colSums(fw_wide), "/") * 100
  fw_rel[is.nan(fw_rel)] <- 0
  fw_plot <- log10(fw_rel + 0.01)
  
  ## 6. annotations ---------------------------------------------------------
  ann_col <- meta %>% filter(sample %in% colnames(fw_plot)) %>%
    column_to_rownames("sample") %>% select(species, pond)
  ann_row <- fw %>% distinct(taxon, group) %>%
    column_to_rownames("taxon")
  ann_row <- ann_row[rownames(fw_plot), , drop = FALSE]
  
  ann_colours <- list(species = species_cols,
                      pond    = pond_cols,
                      group   = group_cols)
  
  ## 7. draw to console (no filename = shows in Plots pane) ------------------
  ph <- pheatmap(
    fw_plot,
    annotation_col    = ann_col,
    annotation_row    = ann_row,
    annotation_colors = ann_colours,
    clustering_method = "average",
    cluster_rows      = nrow(fw_plot) > 2,
    cluster_cols      = ncol(fw_plot) > 2,
    color             = colorRampPalette(brewer.pal(9,"YlGnBu"))(100),
    border_color      = "grey90",
    fontsize_row      = ifelse(LEVEL == "s", 6, 9),
    main = paste0("Freshwater crustacean prey — ", rank_label,
                  " level (log10 % within Crustacea)")
  )
  
  ## 8. Save to file ---------------------------------------------------
  outfile <- paste0("heatmap_freshwater_crustacea_", rank_label, ".png")
  png(outfile, width = 9, height = 6, units = "in", res = 300)
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  dev.off()
  cat("saved ->", outfile, "\n")
  
  invisible(ph)
}

## ===========================================================================
## RUN: a) family b) genus
## ===========================================================================
make_crust_heatmap("f")   # a) FAMILY level
make_crust_heatmap("g")   # b) GENUS level


## ===========================================================================
## EXPORT: numbers behind the heatmaps
## ===========================================================================
export_crust_table <- function(LEVEL = "g") {
  
  rank_label <- c(f = "family", g = "genus")[LEVEL]
  prefix     <- paste0(LEVEL, "__")
  
  tax_long <- all_long %>%
    mutate(last_elem = sub(".*\\|", "", lineage)) %>%
    filter(startsWith(last_elem, prefix)) %>%
    mutate(taxon  = sub(paste0("^", prefix), "", last_elem),
           family = ifelse(grepl("\\|f__", lineage),
                           sub(".*\\|f__([^|]+).*", "\\1", lineage), NA)) %>%
    select(sample, family, taxon, count) %>%
    inner_join(whitelist, by = "family")
  
  ## within-sample % of crustacean reads
  per_sample <- tax_long %>%
    group_by(sample, taxon, group) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(pct = 100 * count / sum(count)) %>%
    ungroup() %>%
    left_join(meta, by = "sample")
  
  ## 1. full per-sample matrix (taxon x sample, % values)
  wide <- per_sample %>%
    select(taxon, sample, pct) %>%
    pivot_wider(names_from = sample, values_from = pct, values_fill = 0) %>%
    arrange(desc(rowSums(across(where(is.numeric)))))
  write_csv(wide, paste0("heatmap_values_", rank_label, "_persample.csv"))
  
  ## 2. mean % by species
  by_species <- per_sample %>%
    group_by(taxon, group, species) %>%
    summarise(mean_pct = mean(pct), .groups = "drop") %>%
    pivot_wider(names_from = species, values_from = mean_pct,
                values_fill = 0, names_prefix = "mean_") %>%
    mutate(diff_brem_minus_aus = mean_Ubrem - mean_Uaus) %>%
    arrange(desc(abs(diff_brem_minus_aus)))
  write_csv(by_species, paste0("heatmap_summary_", rank_label, "_by_species.csv"))
  
  ## 3. mean % by situation
  by_situation <- per_sample %>%
    filter(situation != "other") %>%
    group_by(taxon, group, situation) %>%
    summarise(mean_pct = mean(pct), .groups = "drop") %>%
    pivot_wider(names_from = situation, values_from = mean_pct, values_fill = 0)
  write_csv(by_situation, paste0("heatmap_summary_", rank_label, "_by_situation.csv"))
  
  cat("wrote 2 CSVs for", rank_label, "\n")
  invisible(list(persample = wide, species = by_species, situation = by_situation))
}

export_crust_table("f")
export_crust_table("g")
