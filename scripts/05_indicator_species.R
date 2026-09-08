# Identify ASVs characteristic of each treatment, and draw Figure 4
#  requires data/processed/James_MU42022_filtered_not_normalized.RDS
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)
#
#  IndVal is run twice: pooled across the four matched days for Table S2 and
#  Figure 4, and separately within each day to show when indicators emerged.

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
#install.packages("indicspecies")
library(indicspecies)
#install.packages("ggrepel")
library(ggrepel)

source(here::here("scripts/00_colors.R"))

# User set variables
input.FN     <- "data/processed/James_MU42022_filtered_not_normalized.RDS"
table_s2.FN  <- "tables/TableS2_indval_all_90.csv"
table2.FN    <- "tables/Table2_top5_indicators.csv"
full.FN      <- "tables/indval_matched_days_full.csv"
per_day.FN   <- "tables/indval_per_day_matched_sig.csv"
fig4.FN      <- "figures/Figure4_indval_bubble"
prevalence   <- 0.10
n_perm       <- 9999
alpha        <- 0.05
seed         <- 100

matched_days <- c("Day 01", "Day 03", "Day 06", "Day 15")

# Manuscript value. A mismatch means the input or the filtering has changed.
expected_indicators <- 90L

set.seed(seed)


#### 01. Load and filter ####
stopifnot(file.exists(input.FN))
ps.obj <- readRDS(input.FN)
ps.obj <- subset_samples(ps.obj, !Organism %in% c("Algae", "Chlorophyta", "Diatoms"))
ps.obj <- prune_taxa(taxa_sums(ps.obj) > 0, ps.obj)

md.df <- data.frame(sample_data(ps.obj), stringsAsFactors = FALSE, check.names = FALSE)
md.df <- md.df[, colnames(md.df) != "Sample.ID"]
md.df$Sample.ID <- rownames(md.df)

md.df <- md.df %>%
  mutate(Replicate = str_extract(Sample.ID, pattern = "(?<=-r)\\d+") %>% replace_na("1") %>% factor()
         , Tank      = str_extract(Sample.ID, pattern = "(?<=T)\\d+(?=-)") %>% factor()
         , Day_num   = as.numeric(str_extract(Age, pattern = "\\d+"))
         ) %>%
  filter(Age %in% matched_days)

sample_data(ps.obj) <- sample_data(md.df)
ps.obj  <- prune_samples(sample_names(ps.obj) %in% rownames(md.df), ps.obj)
ps.obj  <- prune_taxa(taxa_sums(ps.obj) > 0, ps.obj)
ps_filt.obj <- filter_taxa(ps.obj, function(x) sum(x > 0) > prevalence * length(x), prune = TRUE)

nsamples(ps_filt.obj)   # 45
ntaxa(ps_filt.obj)      # 739
table(sample_data(ps_filt.obj)$Treatment, sample_data(ps_filt.obj)$Age)

# tax_ind, not tax_df: script 04 defines a tax_df and both are sourced by run_all
tax_ind.df <- as.data.frame(tax_table(ps_filt.obj))

comm.mat <- as.data.frame(otu_table(ps_filt.obj))
if(taxa_are_rows(ps_filt.obj)){ comm.mat <- as.data.frame(t(comm.mat)) }


#### 02. IndVal pooled across matched days ####
treatment_groups <- as.character(sample_data(ps_filt.obj)$Treatment)

# Do not re-seed here or inside the per-day loop below. Every multipatt call
# draws from one stream set in the front matter, and the published counts (90
# pooled, 59 at Day 15) depend on that order. Re-seeding per call shifts ASVs
# sitting near p = 0.05 and changes the Day 15 count to 63.
indval.fit <- multipatt(x = comm.mat, cluster = treatment_groups
                        , func = "r.g", control = how(nperm = n_perm), duleg = TRUE)

# multipatt indexes groups by the sorted order of the cluster vector, so the
# index-to-name map is only valid if that order is what we expect
sorted_levels <- sort(unique(treatment_groups))
stopifnot(identical(sorted_levels
                    , c("Antibiotics", "Antibiotics + HT", "Control", "High temperature")))
idx_map <- c("1" = "Antibiotics", "2" = "Antibiotics + HT"
             , "3" = "Control",     "4" = "High temperature")

add_taxonomy <- function(sign.df){
  sign.df %>%
    rownames_to_column(var = "ASV") %>%
    left_join(tax_ind.df %>%
                rownames_to_column(var = "ASV") %>%
                dplyr::select(ASV, Phylum, Class, Order, Family, Genus)
              , by = "ASV") %>%
    mutate(across(where(is.character), ~replace_na(.x, "Unclassified")))
}

indval_sig.df <- indval.fit$sign %>%
  filter(p.value < alpha) %>%
  add_taxonomy() %>%
  mutate(Treatment = idx_map[as.character(index)]) %>%
  arrange(p.value)

indval_full.df <- indval.fit$sign %>%
  add_taxonomy() %>%
  mutate(Treatment = idx_map[as.character(index)]) %>%
  arrange(p.value)

nrow(indval_sig.df)               # 90
table(indval_sig.df$Treatment)    # 59 AB, 22 AB+HT, 1 Control, 8 HT

if(nrow(indval_sig.df) != expected_indicators){
  warning(sprintf("IndVal count is %d, manuscript reports %d"
                  , nrow(indval_sig.df), expected_indicators))
}


#### 03. Table 2, top five per treatment ####
table2.df <- indval_sig.df %>%
  group_by(Treatment) %>%
  arrange(desc(stat)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  dplyr::select(Taxon = Genus, Family, Class, Treatment, IndVal = stat, p_value = p.value) %>%
  mutate(Taxon   = ifelse(Taxon == "Unclassified", paste0(Family, " (unclassified)"), Taxon)
         , IndVal  = round(IndVal, 3)
         , p_value = round(p_value, 4)
         ) %>%
  arrange(Treatment, desc(IndVal))

as.data.frame(table2.df)


#### 04. IndVal within each day ####
# Shows when indicators emerged. A day needs at least two groups with two or
# more samples for multipatt to run.
per_day.list <- list()

for(day in matched_days){
  ps_day.obj <- prune_samples(sample_data(ps_filt.obj)$Age == day, ps_filt.obj)
  ps_day.obj <- prune_taxa(taxa_sums(ps_day.obj) > 0, ps_day.obj)

  trt_n        <- table(sample_data(ps_day.obj)$Treatment)
  valid_groups <- names(trt_n[trt_n >= 2])
  if(length(valid_groups) < 2){ next }

  ps_day.obj <- prune_samples(sample_data(ps_day.obj)$Treatment %in% valid_groups, ps_day.obj)
  ps_day.obj <- prune_taxa(taxa_sums(ps_day.obj) > 0, ps_day.obj)

  comm_day.mat <- as.data.frame(otu_table(ps_day.obj))
  if(taxa_are_rows(ps_day.obj)){ comm_day.mat <- as.data.frame(t(comm_day.mat)) }

  fit_day <- multipatt(x = comm_day.mat
                       , cluster = as.character(sample_data(ps_day.obj)$Treatment)
                       , func = "r.g", control = how(nperm = n_perm), duleg = TRUE)

  per_day.list[[day]] <- fit_day$sign %>%
    filter(p.value < alpha) %>%
    add_taxonomy() %>%
    mutate(Day = day) %>%
    relocate(Day, .after = p.value) %>%
    arrange(p.value)
}

per_day.df <- bind_rows(per_day.list)
sum(per_day.df$Day == "Day 15")   # 59, the persistence result


#### 05. Figure 4 ####
# Label with Genus where classified, else "Family sp.". A numeric suffix is
# added only when the same label repeats within a treatment.
plot.df <- indval_sig.df %>%
  mutate(Treatment   = factor(Treatment, levels = TREATMENT_ORDER)
         , Taxon_label = if_else(Genus != "Unclassified", Genus, paste0(Family, " sp."))
         ) %>%
  group_by(Treatment, Taxon_label) %>%
  mutate(n_dup = n(), row_n = row_number()
         , Taxon_label = if_else(n_dup > 1, paste0(Taxon_label, " (", row_n, ")"), Taxon_label)
         ) %>%
  ungroup() %>%
  dplyr::select(-n_dup, -row_n)

top5.df <- plot.df %>%
  group_by(Treatment) %>%
  arrange(desc(stat)) %>%
  slice_head(n = 5) %>%
  ungroup()

# Colour is the only encoding for class here, so the number shown is capped at
# what CLASS_PALETTE can supply. Interpolating a qualitative palette past its
# own length produces neighbouring colours that are not separable.
top_cls <- plot.df %>% count(Class, sort = TRUE) %>% slice_head(n = N_CLASSES_SHOWN) %>% pull(Class)
sort(top_cls)
length(setdiff(unique(plot.df$Class), top_cls))   # pooled into Other

# Other is forced last so it takes the grey, not a hue
cls_levels <- c(sort(top_cls), "Other")
plot.df <- plot.df %>% mutate(Class = factor(ifelse(Class %in% top_cls, Class, "Other")
                                             , levels = cls_levels))
top5.df <- top5.df %>% mutate(Class = factor(ifelse(Class %in% top_cls, Class, "Other")
                                             , levels = cls_levels))
cls_colors <- setNames(c(CLASS_PALETTE[seq_along(top_cls)], CLASS_OTHER_COLOR), cls_levels)

fig4.plot <- ggplot(plot.df, aes(x = stat, y = Treatment)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "gray60", linewidth = 0.8) +
  geom_point(aes(size = stat, colour = Class, fill = Class)
             , position = position_jitter(width = 0, height = 0.18, seed = seed)
             , alpha = 0.70, shape = 21, stroke = 0.3
             ) +
  # No position_jitter here. geom_text_repel inherits the already-jittered
  # positions through the data; adding a second jitter misaligns every label.
  geom_text_repel(data = top5.df
                  , aes(label = Taxon_label)
                  , size = 3, fontface = "italic", max.overlaps = 25
                  , box.padding = 0.5, point.padding = 0.3
                  , segment.size = 0.2, segment.color = "gray40"
                  , min.segment.length = 0.1, seed = seed
                  ) +
  scale_colour_manual(values = cls_colors, name = "Taxonomic Class") +
  scale_fill_manual(values   = cls_colors, name = "Taxonomic Class") +
  scale_size_continuous(name = "IndVal", range = c(1, 12)
                        , breaks = c(0.3, 0.4, 0.5, 0.6)
                        , labels = c("0.3", "0.4", "0.5", "0.6")
                        ) +
  # Axis starts just below the smallest significant IndVal so the panel is not
  # mostly empty. Check min(indval_sig.df$stat) if the input data change.
  scale_x_continuous(limits = c(0.30, 0.65), breaks = seq(0.30, 0.65, 0.05), expand = c(0.02, 0)) +
  labs(x = "Indicator Value", y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.line    = element_line(colour = "black", linewidth = 0.8)
        , axis.ticks   = element_line(colour = "black", linewidth = 0.6)
        , axis.text    = element_text(size = 12, colour = "black")
        , axis.title   = element_text(size = 13, face = "bold")
        , legend.position = "right"
        , legend.title    = element_text(size = 11, face = "bold")
        , legend.text     = element_text(size = 10)
        ) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1), ncol = 1)
         , fill   = guide_legend(override.aes = list(size = 5, alpha = 1), ncol = 1)
         , size   = guide_legend(override.aes = list(colour = "black", fill = "gray60"))
         )

min(indval_sig.df$stat)   # 0.362, above the 0.30 axis floor


#### 06. Write outputs ####
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables",  showWarnings = FALSE, recursive = TRUE)

write.csv(indval_sig.df,  file = table_s2.FN, row.names = FALSE)
write.csv(indval_full.df, file = full.FN,     row.names = FALSE)
write.csv(table2.df,      file = table2.FN,   row.names = FALSE)
write.csv(per_day.df,     file = per_day.FN,  row.names = FALSE)

ggsave(filename = paste0(fig4.FN, ".png"), plot = fig4.plot, width = 10, height = 6, dpi = 300)
ggsave(filename = paste0(fig4.FN, ".pdf"), plot = fig4.plot, width = 10, height = 6, device = cairo_pdf)

# Go to scripts/06_composition_barplots.R
