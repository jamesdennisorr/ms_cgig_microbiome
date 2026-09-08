# Summarize class-level community composition, and draw Figure 3
#  requires data/processed/James_MU42022_filtered_not_normalized.RDS
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)
#
#  Relative abundance is taken from raw counts, not a normalized object, because
#  proportions within a sample are what the figure shows.

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

source(here::here("scripts/00_colors.R"))

# User set variables
input.FN <- "data/processed/James_MU42022_filtered_not_normalized.RDS"
fig3.FN  <- "figures/Figure3_class_composition"
seed     <- 100

# All six days are drawn. Control is absent on Days 8 and 10.
keep_days <- c("Day 01", "Day 03", "Day 06", "Day 08", "Day 10", "Day 15")

set.seed(seed)


#### 01. Load and prepare ####
ps.obj <- readRDS(input.FN)
ps.obj <- subset_samples(ps.obj, !Organism %in% c("Algae", "Chlorophyta", "Diatoms"))
ps.obj <- prune_taxa(taxa_sums(ps.obj) > 0, ps.obj)

md.df <- data.frame(sample_data(ps.obj), stringsAsFactors = FALSE, check.names = FALSE)
md.df <- md.df[, colnames(md.df) != "Sample.ID"]
md.df$Sample.ID <- rownames(md.df)

md.df <- md.df %>%
  mutate(Replicate = str_extract(Sample.ID, pattern = "(?<=-r)\\d+") %>% replace_na("1") %>% factor()
         , Tank      = str_extract(Sample.ID, pattern = "(?<=T)\\d+(?=-)") %>% factor()
         ) %>%
  filter(Age %in% keep_days)

sample_data(ps.obj) <- sample_data(md.df)
ps.obj <- prune_samples(sample_names(ps.obj) %in% rownames(md.df), ps.obj)
ps.obj <- prune_taxa(taxa_sums(ps.obj) > 0, ps.obj)

nsamples(ps.obj)   # 61
table(sample_data(ps.obj)$Treatment, sample_data(ps.obj)$Age)


#### 02. Class-level relative abundance ####
ps_rel.obj   <- transform_sample_counts(ps.obj, function(x) x / sum(x))
ps_class.obj <- tax_glom(ps_rel.obj, taxrank = "Class", NArm = FALSE)

class.df <- psmelt(ps_class.obj) %>%
  rename(RelAbundance = Abundance) %>%
  mutate(Class   = ifelse(is.na(Class), paste0("Unclassified_", Phylum), Class)
         , Day_num = as.numeric(str_extract(Age, pattern = "\\d+"))
         )

class_summary.df <- class.df %>%
  group_by(Treatment, Age, Day_num, Class) %>%
  summarise(mean_rel_abund = mean(RelAbundance, na.rm = TRUE)
            , sd_rel_abund   = sd(RelAbundance,   na.rm = TRUE)
            , n              = n()
            , .groups = "drop"
            )

# Values quoted in Results 3.3
class_summary.df %>%
  filter(Class == "Gammaproteobacteria", Age == "Day 01") %>%
  dplyr::select(Treatment, mean_rel_abund, sd_rel_abund)   # AB 0.82, AB+HT 0.73, Ctrl 0.52, HT 0.39


#### 03. Figure 3 ####
# Show the most abundant classes individually, pool the rest. The cap is set by
# CLASS_PALETTE: Wong gives seven hues plus a grey, and asking for more forces
# colours the palette cannot supply.
top_classes <- class_summary.df %>%
  group_by(Class) %>%
  summarise(max_abund = max(mean_rel_abund), .groups = "drop") %>%
  arrange(desc(max_abund)) %>%
  slice_head(n = N_CLASSES_SHOWN) %>%
  pull(Class)

sort(top_classes)
length(setdiff(unique(class_summary.df$Class), top_classes))   # pooled into Other

# Other is forced last so it takes the grey, not a hue
class_levels <- c(sort(top_classes), "Other")
plot.df <- class_summary.df %>%
  mutate(Class_plot = factor(ifelse(Class %in% top_classes, Class, "Other"), levels = class_levels)
         , Treatment  = factor(Treatment, levels = TREATMENT_ORDER)
         ) %>%
  group_by(Treatment, Age, Day_num, Class_plot) %>%
  summarise(mean_rel_abund = sum(mean_rel_abund), .groups = "drop")

class_colors <- setNames(c(CLASS_PALETTE[seq_along(top_classes)], CLASS_OTHER_COLOR), class_levels)

fig3.plot <- ggplot(plot.df, aes(x = factor(Day_num), y = mean_rel_abund, fill = Class_plot)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  facet_wrap(~ Treatment, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Day post-fertilization", y = "Mean relative abundance", fill = "Class") +
  theme_classic(base_size = 12) +
  theme(legend.position  = "right"
        , axis.text.x      = element_text(angle = 0)
        , strip.background = element_rect(fill = "grey95")
        , strip.text       = element_text(face = "bold")
        )


#### 04. Write outputs ####
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

ggsave(filename = paste0(fig3.FN, ".png"), plot = fig3.plot, width = 12, height = 8, dpi = 300)
ggsave(filename = paste0(fig3.FN, ".pdf"), plot = fig3.plot, width = 12, height = 8)

# Go to scripts/07_survival_analysis.R
