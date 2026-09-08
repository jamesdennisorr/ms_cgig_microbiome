# Regenerate every dataset figure quoted in Methods 2.5 and 2.6
#  requires data/processed/James_MU42022_filtered_not_normalized.RDS
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)
#
#  Writes nothing. Exists so that no number in the Methods is untraceable.

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
#install.packages("here")
library(here)
#BiocManager::install("phyloseq")
library(phyloseq)

# User set variables
input.FN     <- "data/processed/James_MU42022_filtered_not_normalized.RDS"
prevalence   <- 0.10     # Methods 2.5, ASV must appear in this fraction of samples
rarefy_depth <- 5444     # Methods 2.6
seed         <- 100

matched_days <- c("Day 01", "Day 03", "Day 06", "Day 15")

set.seed(seed)


#### 01. Apply the exclusions every analysis script uses ####
ps.obj <- readRDS(input.FN)

# Algae are feed, not oyster microbiome. Depths must be measured after this step
# or they describe a dataset no script actually analyses.
ps.obj <- subset_samples(ps.obj, !Organism %in% c("Algae", "Chlorophyta", "Diatoms"))
ps.obj <- prune_taxa(taxa_sums(ps.obj) > 0, ps.obj)

nsamples(ps.obj)   # 61

# Control is absent on Days 8 and 10, so alpha and beta diversity use these four
m.obj <- prune_samples(sample_data(ps.obj)$Age %in% matched_days, ps.obj)
m.obj <- prune_taxa(taxa_sums(m.obj) > 0, m.obj)


#### 02. Methods 2.5 ####
keep   <- rowSums(otu_table(m.obj) > 0) >= prevalence * nsamples(m.obj)
depths <- sample_sums(m.obj)

sum(keep)             # 739 ASVs
nsamples(m.obj)       # 45 samples
min(depths)           # 3909
max(depths)           # 115374
round(median(depths)) # 31447


#### 03. Methods 2.6 ####
# Rarefying below the minimum drops that sample, which is why alpha diversity
# has one fewer sample than beta diversity
sum(depths <  rarefy_depth)   # 1 sample dropped
sum(depths >= rarefy_depth)   # 44, n for alpha diversity
nsamples(m.obj)               # 45, n for beta diversity (CSS, not rarefied)

# End
