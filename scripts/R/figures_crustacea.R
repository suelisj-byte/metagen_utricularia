###############################################################################
##  Main crustacean analysis
##  U. bremii vs U. australis (Aisch region)
###############################################################################

## 0. PACKAGES ================================================================
packages <- c("tidyverse", "vegan", "ggpubr", "patchwork", "FSA", "scales", "ggtext")
to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
lapply(packages, library, character.only = TRUE)

## 1. SETTINGS ================================================================
data_dir  <- "/path/to/mpa_data"
MIN_READS <- 100
set.seed(42)

situation_cols <- c(
  "0_cooccur"   = "#78bc21",  # situation 0 – co-occurrence  (chartreuse)
  "1_bremii"    = "#1976D2",  # situation 1 – U. bremii      (blue)
  "2_australis" = "#b6520f"   # situation 2 – U. australis   (brown orange)
)

species_cols <- c(Ubrem = "#7A3E48", Uaus = "#E6B800")
sp_lab_md <- c(Ubrem = "*U.bremii*", Uaus = "*U. australis*")

pond_cols <- c(
  "I"   = "#78bc21",  # sit 0  (chartreuse)
  "II"  = "#1976D2",  # sit 1  (blue)
  "IV"  = "#64B5F6",  # sit 1  (light blue)
  "III" = "#D16002",  # sit 2  (marmelade)
  "V"   = "#E65100",  # sit 2  (dark orange)
  "VI"  = "#FF6E00",  # sit 2  (hotorange)
  "VII" = "#FAB972"   # sit 2  (calm orange)
)

## 2. WHITELIST =====================
whitelist <- tribble(
  ~family,            ~group,
  "Daphniidae","Cladocera","Chydoridae","Cladocera","Bosminidae","Cladocera",
  "Sididae","Cladocera","Macrothricidae","Cladocera",
  "Cyclopidae","Copepoda","Cyclopettidae","Copepoda","Diaptomidae","Copepoda",
  "Temoridae","Copepoda","Canthocamptidae","Copepoda",
  "Cyprididae","Ostracoda","Candonidae","Ostracoda","Darwinulidae","Ostracoda",
  "Gammaridae","Amphipoda","Crangonyctidae","Amphipoda","Hyalellidae","Amphipoda"
)

## 3. READ .mpa FILES ========================================================
files <- list.files(data_dir, pattern = "^S\\d{2}\\.mpa\\.txt$", full.names = TRUE)
names(files) <- sub("\\.mpa\\.txt$", "", basename(files))
stopifnot(length(files) > 0)

all_long <- imap_dfr(files, ~ read_tsv(.x, col_names = c("lineage","count"),
                                       show_col_types = FALSE) %>%
                       mutate(sample = .y))

## 4. BUILD CRUSTACEAN GENUS MATRIX ==========================================
gen_long <- all_long %>%
  mutate(last_elem = sub(".*\\|", "", lineage)) %>%
  filter(startsWith(last_elem, "g__")) %>%
  mutate(genus  = sub("^g__", "", last_elem),
         family = sub(".*\\|f__([^|]+).*", "\\1", lineage)) %>%
  select(sample, family, genus, count) %>%
  inner_join(whitelist, by = "family")

abund <- gen_long %>%
  group_by(sample, genus) %>% summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = count, values_fill = 0) %>%
  column_to_rownames("sample") %>% as.matrix()
abund <- abund[order(rownames(abund)), ]

## 5. METADATA (Pond VIII pooled into II; III & IV -> "other") ===============
meta <- tibble::tribble(
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

## 6. QUALITY FILTER =========================================================
depth <- rowSums(abund)
keep  <- names(depth)[depth >= MIN_READS]
cat("Excluded (<", MIN_READS, "crustacean reads):",
    paste(setdiff(rownames(abund), keep), collapse = ", "), "\n")

abund <- abund[keep, , drop = FALSE]
abund <- abund[, colSums(abund) > 0, drop = FALSE]
meta  <- meta %>% filter(sample %in% rownames(abund))
meta  <- meta[match(rownames(abund), meta$sample), ]
stopifnot(all(rownames(abund) == meta$sample))

## 7. NORMALISE ==============================================================
rel <- sweep(abund, 1, rowSums(abund), "/")

## 8. ALPHA DIVERSITY (Shannon) ==============================================
alpha <- tibble(
  sample  = rownames(abund),
  Shannon = diversity(abund, index = "shannon"),
  depth   = rowSums(abund)
) %>% left_join(meta, by = "sample")

write_csv(alpha, "alpha_shannon_crustacea.csv")

cat("\n=== Shapiro-Wilk (Shannon) ===\n");            print(shapiro.test(alpha$Shannon))
cat("\n=== Wilcoxon: Shannon ~ species ===\n");        print(wilcox.test(Shannon ~ species, data = alpha))

alpha_sit <- alpha %>% filter(situation != "other") %>% droplevels()
cat("\n=== Kruskal-Wallis: Shannon ~ situation ===\n"); print(kruskal.test(Shannon ~ situation, data = alpha_sit))
# print(FSA::dunnTest(Shannon ~ situation, data = alpha_sit, method = "bh"))  # if sig.

## 8a. FIGURE 1A --------------------------------------------------------------
sit_labels <- c("0_cooccur"   = "Co-occurrence",
                "1_bremii"    = "*U. bremii*",
                "2_australis" = "*U. australis*")

pA1 <- ggplot(alpha, aes(species, Shannon, fill = species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 2) +
  scale_fill_manual(values = species_cols) +
  scale_x_discrete(labels = sp_lab_md) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.4) +
  labs(title = "a) Species", x = NULL, y = "Shannon (H')") +
  guides(fill = "none") + 
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_markdown())

## 8b. Figure 1B: by situation (colour = situation, shape = species)
pA2 <- ggplot(alpha_sit, aes(situation, Shannon, fill = situation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(aes(color = species), width = 0.15, size = 2) +
  scale_fill_manual(values = situation_cols, guide = "none") +
  scale_color_manual(values = species_cols,
                     labels = sp_lab_md,
                     name = "Species") +
  scale_x_discrete(labels = sit_labels) +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  labs(title = "b) Situation", x = "Situation", y = "Shannon (H')") +
  theme_bw(base_size = 12) +
  theme(axis.text.x=element_markdown(), legent.text=element_markdown)

alpha_fig <- pA1 + pA2 + plot_layout(widths = c(1, 1.4))
ggsave("fig1_alpha_diversity_crustacea.png", alpha_fig, width = 11, height = 5, dpi = 300)
print(alpha_fig)

## 9. BETA DIVERSITY (Bray-Curtis) ===========================================
bc <- vegdist(rel, method = "bray")

## 9a. FIGURE 4 — NMDS -------------------------------------------------------
nmds <- metaMDS(bc, k = 2, trymax = 100, trace = FALSE)
cat("\n=== NMDS stress ===\n"); print(nmds$stress)

nmds_df <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("sample") %>% left_join(meta, by = "sample")

pB1 <- ggplot(nmds_df, aes(NMDS1, NMDS2)) +
  stat_ellipse(aes(colour = species, group = species),
               linewidth = 0.6, show.legend = FALSE) +
  geom_point(aes(colour = species, shape = pond), size = 3.2) +
  scale_colour_manual(values = species_cols,
                      labels = sp_lab_md,
                      name = "Species") +
  scale_shape_manual(values = 1:nlevels(meta$pond), name = "Pond", drop = FALSE) +
  labs(title = "NMDS (Bray-Curtis) of crustacean prey") +
  theme_bw(base_size = 12) +
  theme(legend.text = element_markdown())

ggsave("fig4_beta_NMDS_crustacea.png", pB1, width = 8, height = 6, dpi = 300)
print(pB1)

## 9b. PCoA (appendix) — stray-paren bug fixed -------------------------------
pcoa <- cmdscale(bc, k = 2, eig = TRUE)
var_expl <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)
pcoa_df <- as_tibble(pcoa$points, .name_repair = ~c("PCoA1","PCoA2")) %>%
  mutate(sample = rownames(pcoa$points)) %>% left_join(meta, by = "sample")

pB2 <- ggplot(pcoa_df, aes(PCoA1, PCoA2, colour = species, shape = pond)) +
  geom_point(size = 3.2) +
  scale_colour_manual(values = species_cols,                       # <-- fixed
                      labels = sp_lab_md,
                      name = "Species") +
  scale_shape_manual(values = 1:nlevels(meta$pond), name = "Pond", drop = FALSE) +
  labs(title = "PCoA (Bray-Curtis) of crustacean prey",
       x = paste0("PCoA1 (", var_expl[1], "%)"),
       y = paste0("PCoA2 (", var_expl[2], "%)")) +
  theme_bw(base_size = 12) +
  theme(legend.text = element_markdown())

ggsave("appB_beta_PCoA_crustacea.png", pB2, width = 8, height = 6, dpi = 300)
print(pB2)

## 10. PERMANOVA + PERMDISP ==================================================
cat("\n=== PERMANOVA: composition ~ species ===\n")
print(adonis2(bc ~ species, data = meta, permutations = 999))

meta_sit <- meta %>% filter(situation != "other") %>% droplevels()
rel_sit  <- rel[meta_sit$sample, ]
bc_sit   <- vegdist(rel_sit, method = "bray")

cat("\n=== PERMANOVA: composition ~ situation ===\n")
print(adonis2(bc_sit ~ situation, data = meta_sit, permutations = 999))

cat("\n=== PERMDISP: dispersion ~ species ===\n")
print(permutest(betadisper(bc, meta$species), permutations = 999))

cat("\n=== PERMDISP: dispersion ~ situation ===\n")
print(permutest(betadisper(bc_sit, meta_sit$situation), permutations = 999))

## 11. FIGURE 4: relative abundance by genus =====

sample_order <- meta %>% arrange(species, pond) %>% pull(sample)
 
genus_long <- gen_long %>%
  group_by(sample, genus) %>% summarise(count = sum(count), .groups = "drop") %>%
  filter(sample %in% rownames(abund)) %>%          # QC-passing samples only
  left_join(meta, by = "sample") %>%
  mutate(sample = factor(sample, levels = sample_order))
 
ann_df <- meta %>%
  filter(sample %in% sample_order) %>%
  mutate(sample = factor(sample, levels = sample_order))
 
strip_theme <- theme_bw(base_size = 11) +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 9),
        panel.grid   = element_blank())
 
### species colour strip
p_sp <- ggplot(ann_df, aes(sample, 1, fill = species)) +
  geom_tile(width=0.9) +
  scale_fill_manual(values = species_cols, labels = sp_lab_md, name = "Species") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Species") +
  strip_theme + theme(legend.text = element_markdown())
 
### pond colour strip
p_pond <- ggplot(ann_df, aes(sample, 1, fill = pond)) +
  geom_tile(width=0.9) +
  scale_fill_manual(values = pond_cols, name = "Pond", drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Pond") +
  strip_theme
 
## main stacked bars (genus)
p_bars <- ggplot(genus_long, aes(sample, count, fill = genus)) +
  geom_col(position = "fill", width = 0.9) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Sample", y = "Relative read abundance", fill = "Genus") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        panel.grid  = element_blank())
 
## stack strips above the bars; patchwork aligns the x axes and collects legends
fig5 <- p_sp / p_pond / p_bars +
  plot_layout(heights = c(0.5, 0.5, 10), guides = "collect")
 
ggsave("fig4_genus_relabund_crustacea.png", fig5, width = 13, height = 6, dpi = 300)
print(fig5)

## 12. SESSION INFO ==========================================================
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("\nDone. Figures + CSV + sessionInfo written to working directory.\n")
