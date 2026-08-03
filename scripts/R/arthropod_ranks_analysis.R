###############################################################################
##  Utricularia arthropod prey
##  Alpha + beta diversity across taxonomic ranks
##  U. bremii vs U. australis (Aisch region)
##
##  Runs the full analysis at CLASS, ORDER and GENUS level in one pass and
##  collects every test statistic into a single tidy table.
###############################################################################

## 0. PACKAGES ================================================================
packages <- c("tidyverse", "vegan", "ggpubr", "patchwork", "FSA", "ggtext")
to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
lapply(packages, library, character.only = TRUE)

## 1. SETTINGS ================================================================
data_dir  <- "/path/to/mpa_data"
MIN_READS <- 100
RANKS     <- c("c", "o", "f", "g")      # class, order, family, genus
set.seed(42)

rank_names <- c(c = "Class", o = "Order", f = "Family", g = "Genus")

## --- palettes (identical to the crustacean script) --------------------------
situation_cols <- c(
  "0_cooccur"   = "#78bc21",   # situation 0 – co-occurrence
  "1_bremii"    = "#1976D2",   # situation 1 – U. bremii
  "2_australis" = "#E65100"    # situation 2 – U. australis
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

species_cols <- c(Ubrem = "#7A3E48", Uaus = "#E6B800")

sp_lab_md  <- c(Ubrem = "*U. bremii*", Uaus = "*U. australis*")
sit_labels <- c("0_cooccur"   = "Co-occurrence",
                "1_bremii"    = "*U. bremii*",
                "2_australis" = "*U. australis*")

## 2. READ .mpa FILES (once, reused for every rank) ==========================
files <- list.files(data_dir, pattern = "^S\\d{2}\\.mpa\\.txt$", full.names = TRUE)
names(files) <- sub("\\.mpa\\.txt$", "", basename(files))
stopifnot(length(files) > 0)

all_long <- imap_dfr(files, ~ read_tsv(.x, col_names = c("lineage","count"),
                                       show_col_types = FALSE) %>%
                       mutate(sample = .y))

## 3. METADATA (Pond VIII pooled into II; III & IV -> "other") ===============
meta_master <- tibble::tribble(
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
    pond      = factor(pond, levels = c("I","II","III","IV","V","VI","VII")),
    species   = factor(species, levels = c("Ubrem","Uaus")),
    situation = case_when(
      pond == "I"                 ~ "0_cooccur",
      pond == "II"                ~ "1_bremii",
      pond %in% c("V","VI","VII") ~ "2_australis",
      TRUE                        ~ "other"
    ) %>% factor(levels = c("0_cooccur","1_bremii","2_australis","other"))
  )

## ===========================================================================
## Run the whole analysis at one rank
## ===========================================================================
run_rank <- function(RANK) {
  
  rank_label <- rank_names[RANK]
  prefix     <- paste0(RANK, "__")
  cat("\n\n###########################################################\n")
  cat("##  RANK:", rank_label, "\n")
  cat("###########################################################\n")
  
  ## --- abundance matrix ----------------------------------------------------
  abund <- all_long %>%
    mutate(last_elem = sub(".*\\|", "", lineage)) %>%
    filter(startsWith(last_elem, prefix)) %>%
    mutate(taxon = sub(paste0("^", prefix), "", last_elem)) %>%
    group_by(sample, taxon) %>% summarise(count = sum(count), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = count, values_fill = 0) %>%
    column_to_rownames("sample") %>% as.matrix()
  abund <- abund[order(rownames(abund)), , drop = FALSE]
  
  cat("Samples:", nrow(abund), "| Taxa:", ncol(abund), "\n")
  cat("Read range:", min(rowSums(abund)), "-", max(rowSums(abund)), "\n")
  
  ## --- quality filter ------------------------------------------------------
  depth <- rowSums(abund)
  keep  <- names(depth)[depth >= MIN_READS]
  cat("Excluded (<", MIN_READS, "reads):",
      paste(setdiff(rownames(abund), keep), collapse = ", "), "\n")
  
  abund <- abund[keep, , drop = FALSE]
  abund <- abund[, colSums(abund) > 0, drop = FALSE]
  meta  <- meta_master %>% filter(sample %in% rownames(abund))
  meta  <- meta[match(rownames(abund), meta$sample), ]
  stopifnot(all(rownames(abund) == meta$sample))
  
  if (ncol(abund) < 3)
    warning("Only ", ncol(abund), " taxa at ", rank_label,
            " level — diversity metrics will be very coarse.")
  
  ## --- normalise -----------------------------------------------------------
  rel <- sweep(abund, 1, rowSums(abund), "/")
  
  ## --- ALPHA ---------------------------------------------------------------
  alpha <- tibble(sample  = rownames(abund),
                  Shannon = diversity(abund, index = "shannon"),
                  depth   = rowSums(abund)) %>%
    left_join(meta, by = "sample")
  
  write_csv(alpha, paste0("alpha_shannon_", rank_label, ".csv"))
  
  sw <- shapiro.test(alpha$Shannon)
  cat("\n=== Shapiro-Wilk (Shannon) ===\n"); print(sw)
  
  wt <- wilcox.test(Shannon ~ species, data = alpha)
  cat("\n=== Wilcoxon: Shannon ~ species ===\n"); print(wt)
  
  alpha_sit <- alpha %>% filter(situation != "other") %>% droplevels()
  kw <- kruskal.test(Shannon ~ situation, data = alpha_sit)
  cat("\n=== Kruskal-Wallis: Shannon ~ situation ===\n"); print(kw)
  if (kw$p.value < 0.05)
    print(FSA::dunnTest(Shannon ~ situation, data = alpha_sit, method = "bh"))
  
  ## --- FIGURE A1 ------------------------------------------------------------
  pA1 <- ggplot(alpha, aes(species, Shannon, fill = species)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 2) +
    scale_fill_manual(values = species_cols) +
    scale_x_discrete(labels = sp_lab_md) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.4) +
    labs(title = "A) Shannon by species", x = NULL, y = "Shannon (H')") +
    guides(fill = "none") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_markdown())
  
  pA2 <- ggplot(alpha_sit, aes(situation, Shannon, fill = situation)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(aes(colour = species), width = 0.15, size = 2) +
    scale_fill_manual(values = situation_cols, guide = "none") +
    scale_colour_manual(values = species_cols, labels = sp_lab_md, name = "Species") +
    scale_x_discrete(labels = sit_labels) +
    stat_compare_means(method = "kruskal.test", label = "p.format") +
    labs(title = "B) Shannon by situation", x = "Situation", y = "Shannon (H')") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_markdown(), legend.text = element_markdown())
  
  alpha_fig <- pA1 + pA2 + plot_layout(widths = c(1, 1.4))
  ggsave(paste0("alpha_diversity_a1_", rank_label, ".png"),
         alpha_fig, width = 11, height = 5, dpi = 300)
  print(alpha_fig)
  
  ## --- FIGURE B2a ----------------------------------------------------------------
  bc <- vegdist(rel, method = "bray")
  
  nmds <- metaMDS(bc, k = 2, trymax = 100, trace = FALSE)
  cat("\n=== NMDS stress ===\n"); print(nmds$stress)
  
  nmds_df <- as.data.frame(scores(nmds, display = "sites")) %>%
    rownames_to_column("sample") %>% left_join(meta, by = "sample")
  
  pB1 <- ggplot(nmds_df, aes(NMDS1, NMDS2)) +
    stat_ellipse(aes(colour = species, group = species),
                 linewidth = 0.6, show.legend = FALSE) +
    geom_point(aes(colour = species, shape = pond), size = 3.2) +
    scale_colour_manual(values = species_cols, labels = sp_lab_md, name = "Species") +
    scale_shape_manual(values = 1:nlevels(meta$pond), name = "Pond", drop = FALSE) +
    labs(title = paste0("NMDS (Bray-Curtis) — ", rank_label, " level"),
         subtitle = paste("Stress =", round(nmds$stress, 3))) +
    theme_bw(base_size = 12) +
    theme(legend.text = element_markdown())
  
  ggsave(paste0("beta_NMDS_b2a_", rank_label, ".png"), pB1, width = 8, height = 6, dpi = 300)
  print(pB1)
  
  pcoa <- cmdscale(bc, k = 2, eig = TRUE)
  var_expl <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)
  pcoa_df <- as_tibble(pcoa$points, .name_repair = ~c("PCoA1","PCoA2")) %>%
    mutate(sample = rownames(pcoa$points)) %>% left_join(meta, by = "sample")
  
  pB2 <- ggplot(pcoa_df, aes(PCoA1, PCoA2, colour = species, shape = pond)) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = species_cols, labels = sp_lab_md, name = "Species") +
    scale_shape_manual(values = 1:nlevels(meta$pond), name = "Pond", drop = FALSE) +
    labs(title = paste0("PCoA (Bray-Curtis) — ", rank_label, " level"),
         x = paste0("PCoA1 (", var_expl[1], "%)"),
         y = paste0("PCoA2 (", var_expl[2], "%)")) +
    theme_bw(base_size = 12) +
    theme(legend.text = element_markdown())
  
  ggsave(paste0("beta_PCoA_b2b_", rank_label, ".png"), pB2, width = 8, height = 6, dpi = 300)
  print(pB2)
  
  ## --- PERMANOVA + PERMDISP (Table A1-3) ------------------------------------------------
  cat("\n=== PERMANOVA: composition ~ species ===\n")
  pn_sp <- adonis2(bc ~ species, data = meta, permutations = 999); print(pn_sp)
  
  meta_sit <- meta %>% filter(situation != "other") %>% droplevels()
  rel_sit  <- rel[meta_sit$sample, , drop = FALSE]
  bc_sit   <- vegdist(rel_sit, method = "bray")
  
  cat("\n=== PERMANOVA: composition ~ situation ===\n")
  pn_si <- adonis2(bc_sit ~ situation, data = meta_sit, permutations = 999); print(pn_si)
  
  cat("\n=== PERMDISP: dispersion ~ species ===\n")
  pd_sp <- permutest(betadisper(bc, meta$species), permutations = 999); print(pd_sp)
  
  cat("\n=== PERMDISP: dispersion ~ situation ===\n")
  pd_si <- permutest(betadisper(bc_sit, meta_sit$situation), permutations = 999); print(pd_si)
  
  ## --- POND I (unconfounded contrast) --------------------------------------
  cat("\n########## POND I (co-occurrence) only ##########\n")
  pondI <- meta %>% filter(pond == "I")
  pI_w_stat <- pI_w_p <- pI_F <- pI_R2 <- pI_p <- NA_real_
  
  if (nlevels(droplevels(pondI$species)) == 2 && all(table(pondI$species) >= 2)) {
    n_tab <- table(droplevels(pondI$species))
    cat("Pond I n:", paste(names(n_tab), n_tab, collapse = " / "), "\n")
    cat("NOTE: minimum achievable two-sided Wilcoxon p =",
        round(2 / choose(sum(n_tab), min(n_tab)), 4), "\n")
    
    wI <- wilcox.test(Shannon ~ species, data = alpha %>% filter(pond == "I"))
    cat("\n=== Pond I: Shannon ~ species (Wilcoxon) ===\n"); print(wI)
    pI_w_stat <- unname(wI$statistic); pI_w_p <- wI$p.value
    
    relI <- rel[pondI$sample, , drop = FALSE]
    relI <- relI[, colSums(relI) > 0, drop = FALSE]
    bcI  <- vegdist(relI, "bray")
    cat("\n=== Pond I: composition ~ species (PERMANOVA) ===\n")
    aI <- adonis2(bcI ~ species, data = pondI, permutations = 999); print(aI)
    pI_F <- aI$F[1]; pI_R2 <- aI$R2[1]; pI_p <- aI$`Pr(>F)`[1]
  } else {
    cat("Not enough replication in Pond I for a within-pond test.\n")
  }
  
  ## --- collect stats -------------------------------------------------------
  tibble(
    rank              = rank_label,
    n_samples         = nrow(abund),
    n_taxa            = ncol(abund),
    shapiro_p         = sw$p.value,
    wilcox_W          = unname(wt$statistic),
    wilcox_p          = wt$p.value,
    kruskal_chisq     = unname(kw$statistic),
    kruskal_df        = unname(kw$parameter),
    kruskal_p         = kw$p.value,
    permanova_sp_F    = pn_sp$F[1],
    permanova_sp_R2   = pn_sp$R2[1],
    permanova_sp_p    = pn_sp$`Pr(>F)`[1],
    permanova_sit_F   = pn_si$F[1],
    permanova_sit_R2  = pn_si$R2[1],
    permanova_sit_p   = pn_si$`Pr(>F)`[1],
    permdisp_sp_F     = pd_sp$tab$F[1],
    permdisp_sp_p     = pd_sp$tab$`Pr(>F)`[1],
    permdisp_sit_F    = pd_si$tab$F[1],
    permdisp_sit_p    = pd_si$tab$`Pr(>F)`[1],
    nmds_stress       = nmds$stress,
    pondI_wilcox_W    = pI_w_stat,
    pondI_wilcox_p    = pI_w_p,
    pondI_permanova_F = pI_F,
    pondI_permanova_R2 = pI_R2,
    pondI_permanova_p = pI_p
  )
}

## ===========================================================================
## RUN ALL RANKS + WRITE COMBINED STATISTICS TABLE
## ===========================================================================
stats_all <- map_dfr(RANKS, run_rank)

write_csv(stats_all, "stats_summary_all_ranks.csv")
cat("\n\n=== COMBINED STATISTICS (all ranks) ===\n")
print(as.data.frame(stats_all))

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("\nDone. Figures + per-rank CSVs + stats_summary_all_ranks.csv written.\n")