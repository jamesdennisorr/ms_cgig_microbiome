# Test whether early-life perturbation changed alpha diversity, and draw Figure 1
#  requires data/processed/James_MU42022_filtered_not_normalized.RDS
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
#install.packages("here")
library(here)
#BiocManager::install("phyloseq")
library(phyloseq)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("rstatix")
library(rstatix)
#install.packages("patchwork")
library(patchwork)

source(here::here("scripts/00_colors.R"))

# User set variables
input.FN     <- "data/processed/James_MU42022_filtered_not_normalized.RDS"
fig1.FN      <- "figures/Figure1_alpha_diversity"
table_s1.FN  <- "tables/TableS1_alpha_diversity.csv"
rarefy_depth <- 5444                # minimum retained library size
seed         <- 100

# Days on which all four treatments were sampled. Control is absent on Days 8
# and 10, so tests are restricted to these four; the figure still shows all six.
matched_days <- c("Day 01", "Day 03", "Day 06", "Day 15")

set.seed(seed)


#### 01. Load and rarefy ####
ps.obj <- readRDS(input.FN)

# Algae are feed, not oyster microbiome
ps.obj <- subset_samples(ps.obj, !Organism %in% c("Algae", "Chlorophyta", "Diatoms"))
ps.obj <- prune_taxa(taxa_sums(ps.obj) > 0, ps.obj)

nsamples(ps.obj)          # 61
min(sample_sums(ps.obj))  # 3909

# Rarefying below the minimum drops that sample, which is why alpha diversity
# runs on one fewer sample than beta diversity
ps_rare.obj <- rarefy_even_depth(ps.obj
                                 , sample.size = rarefy_depth
                                 , rngseed     = seed
                                 , replace     = FALSE
                                 )
nsamples(ps_rare.obj)     # 60

md.df <- data.frame(sample_data(ps_rare.obj), check.names = FALSE)
table(md.df$Treatment, md.df$Age)


#### 02. Diversity indices ####
alpha.df <- estimate_richness(ps_rare.obj, measures = c("Observed", "Shannon", "Simpson")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample.ID") %>%
  # phyloseq turns hyphens into periods in row names, so put them back
  mutate(Sample.ID = str_replace_all(Sample.ID, pattern = "\\.", replacement = "-"))

alpha_all.df <- alpha.df %>%
  left_join(md.df, by = "Sample.ID") %>%
  mutate(Day_num   = as.numeric(str_extract(Age, pattern = "\\d+"))
         , Treatment = factor(Treatment, levels = TREATMENT_ORDER)
         )

nrow(alpha_all.df)        # 60
stopifnot(nrow(alpha_all.df) == nsamples(ps_rare.obj))


#### 03. Kruskal-Wallis and Dunn tests ####
alpha_matched.df <- alpha_all.df %>% filter(Age %in% matched_days)
nrow(alpha_matched.df)    # 44

# One helper per index, so adding an index does not mean copying a block
run_kw <- function(data.df, metric){
  data.df %>%
    group_by(Age) %>%
    kruskal_test(as.formula(paste(metric, "~ Treatment"))) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance(p.col = "p.adj") %>%
    ungroup()
}

run_dunn <- function(data.df, metric){
  data.df %>%
    group_by(Age) %>%
    dunn_test(as.formula(paste(metric, "~ Treatment")), p.adjust.method = "bonferroni") %>%
    add_significance(p.col = "p.adj") %>%
    ungroup() %>%
    mutate(Metric = metric)
}

kw_shannon.df   <- run_kw(alpha_matched.df, metric = "Shannon")
kw_simpson.df   <- run_kw(alpha_matched.df, metric = "Simpson")
kw_observed.df  <- run_kw(alpha_matched.df, metric = "Observed")

dunn_shannon.df  <- run_dunn(alpha_matched.df, metric = "Shannon")
dunn_simpson.df  <- run_dunn(alpha_matched.df, metric = "Simpson")
dunn_observed.df <- run_dunn(alpha_matched.df, metric = "Observed")

# No Kruskal-Wallis test is significant after correction
kw_shannon.df
kw_simpson.df
kw_observed.df

# Only two Dunn comparisons survive, both between perturbed groups on Day 15
table_s1.df <- bind_rows(dunn_shannon.df, dunn_simpson.df, dunn_observed.df)
table_s1.df %>% filter(p.adj < 0.05)   # Shannon AB vs HT, Simpson HT vs AB+HT, both p.adj = 0.028


#### 04. Figure 1 ####
# Means across all six days. Days 8 and 10 are drawn but were not tested.
alpha_summary.df <- alpha_all.df %>%
  group_by(Treatment, Day_num) %>%
  summarise(Shannon_mean  = mean(Shannon,  na.rm = TRUE)
            , Shannon_sd    = sd(Shannon,    na.rm = TRUE)
            , Simpson_mean  = mean(Simpson,  na.rm = TRUE)
            , Simpson_sd    = sd(Simpson,    na.rm = TRUE)
            , Observed_mean = mean(Observed, na.rm = TRUE)
            , Observed_sd   = sd(Observed,   na.rm = TRUE)
            , .groups = "drop"
            )

# One panel builder, called three times. The grey band marks Days 8 to 10,
# where Control is missing and no test was run.
alpha_panel <- function(summary.df, mean_col, sd_col, y_lab, panel_title, show_legend){
  ggplot(summary.df, aes(x = Day_num
                         , y        = .data[[mean_col]]
                         , color    = Treatment
                         , fill     = Treatment
                         , shape    = Treatment
                         , linetype = Treatment
                         , group    = Treatment
                         )) +
    annotate(geom = "rect", xmin = 7, xmax = 11, ymin = -Inf, ymax = Inf
             , fill = "grey90", alpha = 0.6) +
    geom_ribbon(aes(ymin = .data[[mean_col]] - .data[[sd_col]]
                    , ymax = .data[[mean_col]] + .data[[sd_col]])
                , alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values    = TREATMENT_COLORS) +
    scale_fill_manual(values     = TREATMENT_COLORS) +
    scale_shape_manual(values    = TREATMENT_SHAPES) +
    scale_linetype_manual(values = TREATMENT_LINETYPES) +
    scale_x_continuous(breaks = c(1, 3, 6, 8, 10, 15)) +
    labs(x = "Day post-fertilization", y = y_lab, title = panel_title
         , color = NULL, fill = NULL, shape = NULL, linetype = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position  = if(show_legend) "top" else "none"
          , legend.direction = "horizontal"
          , legend.title     = element_blank()
          , plot.title       = element_text(face = "bold", hjust = 0.5)
          )
}

fig1.plot <- alpha_panel(alpha_summary.df, "Shannon_mean",  "Shannon_sd"
                         , y_lab = "Shannon diversity", panel_title = "A  Shannon Diversity"
                         , show_legend = TRUE) /
             alpha_panel(alpha_summary.df, "Simpson_mean",  "Simpson_sd"
                         , y_lab = "Simpson diversity", panel_title = "B  Simpson Diversity"
                         , show_legend = FALSE) /
             alpha_panel(alpha_summary.df, "Observed_mean", "Observed_sd"
                         , y_lab = "Observed ASVs", panel_title = "C  Observed Richness"
                         , show_legend = FALSE)


#### 05. Write outputs ####
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables",  showWarnings = FALSE, recursive = TRUE)

write.csv(table_s1.df, file = table_s1.FN, row.names = FALSE)

ggsave(filename = paste0(fig1.FN, ".png"), plot = fig1.plot, width = 8, height = 10, dpi = 300)
ggsave(filename = paste0(fig1.FN, ".pdf"), plot = fig1.plot, width = 8, height = 10)

# Go to scripts/03_beta_diversity.R
