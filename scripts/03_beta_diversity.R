# Test whether community composition differed among treatments, and draw Figure 2
#  requires data/processed/James_MU42022_filtered_not_normalized.RDS
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)
#
#  Restricted to Days 1, 3, 6 and 15, the days on which all four treatments were
#  sampled. Ordination and tests use the same 45 samples so the figure and the
#  statistics describe one dataset.
#
#  Unweighted UniFrac appears in Table S3 only. Its dispersion is unequal between
#  treatments (betadisper F = 4.74, p = 0.0061), so PERMANOVA on it is not
#  interpretable and no figure is drawn.

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
#install.packages("here")
library(here)
#BiocManager::install("phyloseq")
library(phyloseq)
#install.packages("vegan")
library(vegan)
#install.packages("phangorn")
library(phangorn)
#BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
#install.packages("tidyverse")
library(tidyverse)

source(here::here("scripts/00_colors.R"))

# User set variables
input.FN     <- "data/processed/James_MU42022_filtered_not_normalized.RDS"
fig2.FN      <- "figures/Figure2_BrayCurtis_PCoA"
fig_s1.FN    <- "figures/FigureS1_WeightedUniFrac_PCoA"
table_s3.FN  <- "tables/TableS3_pairwise_permanova.csv"
n_perm       <- 9999
seed         <- 100

matched_days <- c("Day 01", "Day 03", "Day 06", "Day 15")

# Point shapes by sampling day
day_shapes <- c("Day 01" = 16   # circle
                , "Day 03" = 17   # triangle
                , "Day 06" = 15   # square
                , "Day 15" = 9    # diamond plus
                )

set.seed(seed)


#### 01. Load and normalize ####
ps.obj <- readRDS(input.FN)
ps.obj <- subset_samples(ps.obj, !Organism %in% c("Algae", "Chlorophyta", "Diatoms"))
ps.obj <- prune_taxa(taxa_sums(ps.obj) > 0, ps.obj)

nsamples(ps.obj)   # 61

# CSS rather than rarefaction, so no sample is discarded for low depth
mr.obj  <- phyloseq_to_metagenomeSeq(ps.obj)
mr.obj  <- cumNorm(mr.obj, p = cumNormStatFast(mr.obj))
css.obj <- ps.obj
otu_table(css.obj) <- otu_table(MRcounts(mr.obj, norm = TRUE, log = FALSE), taxa_are_rows = TRUE)
css.obj <- prune_samples(sample_sums(css.obj) > 0, css.obj)

# Restrict before ordinating, so the figure and the tests see the same samples
css.obj <- prune_samples(sample_data(css.obj)$Age %in% matched_days, css.obj)
css.obj <- prune_taxa(taxa_sums(css.obj) > 0, css.obj)

nsamples(css.obj)  # 45

# The SILVA tree is unrooted. UniFrac needs a root, and phyloseq picks one at
# random when none is given, so distances change on every call. Midpoint rooting
# fixes it deterministically: no outgroup is available for this dataset.
ape::is.rooted(phy_tree(css.obj))                      # FALSE
phy_tree(css.obj) <- phangorn::midpoint(phy_tree(css.obj))
ape::is.rooted(phy_tree(css.obj))                      # TRUE

md.df <- data.frame(sample_data(css.obj), check.names = FALSE)
md.df <- md.df[, colnames(md.df) != "Sample.ID"]   # phyloseq may already carry one
md.df <- rownames_to_column(md.df, var = "Sample.ID")


#### 02. Distances ####
dist_bray  <- phyloseq::distance(css.obj, method = "bray")
dist_wuni  <- UniFrac(css.obj, weighted = TRUE,  normalized = TRUE)
dist_ununi <- UniFrac(css.obj, weighted = FALSE, normalized = TRUE)


#### 03. Ordination ####
# One builder, called for each metric
pcoa_panel <- function(dist.obj, meta.df, plot_title){
  ord.obj <- ordinate(css.obj, method = "PCoA", distance = dist.obj)
  pct     <- round(100 * ord.obj$values$Relative_eig[1:2], 1)

  plot.df <- data.frame(ord.obj$vectors[, 1:2]) %>%
    rownames_to_column(var = "Sample.ID") %>%
    rename(PC1 = Axis.1, PC2 = Axis.2) %>%
    left_join(meta.df, by = "Sample.ID") %>%
    mutate(Treatment = factor(Treatment, levels = TREATMENT_ORDER)
           , Age       = factor(Age, levels = sort(unique(Age)))
           )

  ggplot(plot.df, aes(x = PC1, y = PC2, color = Treatment, shape = Age)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = TREATMENT_COLORS, name = "Treatment") +
    scale_shape_manual(values = day_shapes,       name = "Day") +
    labs(x     = sprintf("PC1 (%.1f%%)", pct[1])
         , y     = sprintf("PC2 (%.1f%%)", pct[2])
         , title = plot_title
         ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "right"
          , plot.title      = element_text(face = "bold")
          )
}

fig2.plot   <- pcoa_panel(dist_bray, md.df, plot_title = "Bray-Curtis PCoA")
fig_s1.plot <- pcoa_panel(dist_wuni, md.df, plot_title = "Weighted UniFrac PCoA")

fig2.plot$labels$x     # PC1 (23.5%)
fig_s1.plot$labels$x   # PC1 (45.4%)


#### 04. PERMANOVA and dispersion ####
# Re-seeded so each permutation block reproduces on its own, whatever ran before
set.seed(seed)

perm_bray <- adonis2(dist_bray ~ Treatment * Age, data = md.df
                     , permutations = n_perm, by = "margin")
perm_wuni <- adonis2(dist_wuni ~ Treatment * Age, data = md.df
                     , permutations = n_perm, by = "margin")

perm_bray   # Treatment:Age  df 9  R2 0.270  F 2.35  p 1e-04
perm_wuni   # Treatment:Age  df 9  R2 0.267  F 2.70  p 2e-04

# Equal dispersion is a precondition for reading PERMANOVA as a location effect
set.seed(seed)
disp_bray  <- permutest(betadisper(dist_bray,  md.df$Treatment), permutations = n_perm)
disp_wuni  <- permutest(betadisper(dist_wuni,  md.df$Treatment), permutations = n_perm)
disp_ununi <- permutest(betadisper(dist_ununi, md.df$Treatment), permutations = n_perm)

disp_bray    # F 0.37  p 0.78
disp_wuni    # F 0.16  p 0.923
disp_ununi   # F 4.74  p 0.0061, which is why unweighted UniFrac gets no figure


#### 05. Table S3, pairwise PERMANOVA ####
run_pairwise <- function(dist.obj, meta.df){
  combos <- combn(unique(meta.df$Treatment), m = 2, simplify = FALSE)
  raw_p  <- sapply(combos, function(pair){
    sub.df <- filter(meta.df, Treatment %in% pair)
    sub_dm <- as.dist(as.matrix(dist.obj)[sub.df$Sample.ID, sub.df$Sample.ID])
    adonis2(sub_dm ~ Treatment, data = sub.df, permutations = n_perm)$`Pr(>F)`[1]
  })
  names(raw_p) <- sapply(combos, paste, collapse = " vs ")
  p.adjust(raw_p, method = "bonferroni")
}

set.seed(seed)
pw_bray  <- run_pairwise(dist_bray,  md.df)
pw_wuni  <- run_pairwise(dist_wuni,  md.df)
pw_ununi <- run_pairwise(dist_ununi, md.df)

round(pw_bray, 4)    # none significant
round(pw_wuni, 4)    # none significant

table_s3.df <- bind_rows(
  data.frame(Metric = "Bray-Curtis",        Comparison = names(pw_bray),  p_adj_bonf = pw_bray,  row.names = NULL)
  , data.frame(Metric = "Weighted UniFrac",   Comparison = names(pw_wuni),  p_adj_bonf = pw_wuni,  row.names = NULL)
  , data.frame(Metric = "Unweighted UniFrac", Comparison = names(pw_ununi), p_adj_bonf = pw_ununi, row.names = NULL)
  )


#### 06. Write outputs ####
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables",  showWarnings = FALSE, recursive = TRUE)

write.csv(table_s3.df, file = table_s3.FN, row.names = FALSE)

ggsave(filename = paste0(fig2.FN, ".png"),   plot = fig2.plot,   width = 8, height = 6, dpi = 300)
ggsave(filename = paste0(fig2.FN, ".pdf"),   plot = fig2.plot,   width = 8, height = 6)
ggsave(filename = paste0(fig_s1.FN, ".png"), plot = fig_s1.plot, width = 8, height = 6, dpi = 300)
ggsave(filename = paste0(fig_s1.FN, ".pdf"), plot = fig_s1.plot, width = 8, height = 6)

# Go to scripts/04_differential_abundance.R
