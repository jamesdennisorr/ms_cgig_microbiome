# Identify differentially abundant ASVs between Control and each treatment
#  requires data/processed/James_MU42022_filtered_not_normalized.RDS
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)
#
#  Two methods are run side by side because they disagree by construction.
#  ALDEx2 tests centred log-ratio abundances with Welch's t; ANCOM-BC2 estimates
#  a bias-corrected log fold change and flags structural zeros. Table 1 reports
#  the union, not the intersection.

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
# ALDEx2 and ANCOMBC pull in MASS, whose select() masks dplyr::select. Load them
# before tidyverse, and still qualify every select() call as dplyr::select.
#install.packages("here")
library(here)
#BiocManager::install("ALDEx2")
library(ALDEx2)
#BiocManager::install("ANCOMBC")
library(ANCOMBC)
#BiocManager::install("phyloseq")
library(phyloseq)
#install.packages("tidyverse")
library(tidyverse)

source(here::here("scripts/00_colors.R"))

# User set variables
input.FN     <- "data/processed/James_MU42022_filtered_not_normalized.RDS"
table1.FN    <- "tables/Table1_differential_abundance.csv"
table_s4.FN  <- "tables/TableS4_DA_complete.csv"
prevalence   <- 0.10     # ASV must appear in at least this fraction of samples
mc_samples   <- 128      # ALDEx2 Monte Carlo instances
alpha        <- 0.05
seed         <- 100

matched_days <- c("Day 01", "Day 03", "Day 06", "Day 15")

set.seed(seed)


#### 01. Load and filter ####
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
  filter(Age %in% matched_days)

sample_data(ps.obj) <- sample_data(md.df)
ps.obj <- prune_samples(sample_names(ps.obj) %in% rownames(md.df), ps.obj)
ps.obj <- prune_taxa(taxa_sums(ps.obj) > 0, ps.obj)
ps.obj <- filter_taxa(ps.obj, function(x) sum(x > 0) > prevalence * length(x), prune = TRUE)

nsamples(ps.obj)   # 45
ntaxa(ps.obj)      # 739
table(sample_data(ps.obj)$Treatment, sample_data(ps.obj)$Age)

# Control first, so every model contrasts against it
sample_data(ps.obj)$Treatment <- factor(sample_data(ps.obj)$Treatment, levels = TREATMENT_ORDER)

tax.df <- as.data.frame(tax_table(ps.obj)) %>%
  rownames_to_column(var = "ASV") %>%
  dplyr::select(ASV, Phylum, Class, Order, Family, Genus)

# One entry per Control-versus-treatment contrast
contrasts.list <- list(list(trt = "Antibiotics",      label = "Antibiotics vs Control")
                       , list(trt = "High temperature", label = "High temperature vs Control")
                       , list(trt = "Antibiotics + HT", label = "Antibiotics + HT vs Control")
                       )


#### 02. ALDEx2 ####
run_aldex_pair <- function(ps.obj, trt, comparison_label){
  # prune_samples, not subset_samples: the latter evaluates its condition in the
  # caller's frame and does not see a local variable named trt
  keep    <- sample_data(ps.obj)$Treatment %in% c("Control", trt)
  ps_sub  <- prune_samples(keep, ps.obj)
  ps_sub  <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

  otu.mat <- as.data.frame(otu_table(ps_sub))
  if(!taxa_are_rows(ps_sub)){ otu.mat <- as.data.frame(t(otu.mat)) }
  conds   <- as.character(sample_data(ps_sub)$Treatment)

  set.seed(seed)
  res.df <- aldex(reads = otu.mat, conditions = conds, mc.samples = mc_samples
                  , denom = "all", test = "t", effect = TRUE)
  res.df <- as.data.frame(res.df) %>% rownames_to_column(var = "ASV")

  # ALDEx2 orders conditions alphabetically, so the sign of diff.btw flips
  # depending on whether the treatment name sorts before or after "Control".
  # Read direction from the rab.win columns instead; they are named by group.
  rab_trt  <- res.df[[paste0("rab.win.", trt)]]
  rab_ctrl <- res.df[["rab.win.Control"]]
  stopifnot(!is.null(rab_trt), !is.null(rab_ctrl))

  res.df %>%
    mutate(clr_trt_minus_ctrl = rab_trt - rab_ctrl   # positive = enriched in treatment
           , Comparison         = comparison_label
           )
}

aldex_all.df <- bind_rows(lapply(contrasts.list
                                 , function(x) run_aldex_pair(ps.obj, x$trt, x$label)))

# No effect-size filter. Adding one drops ASV93 and Table 1 no longer matches.
aldex_sig.df <- aldex_all.df %>% filter(we.eBH < alpha)

nrow(aldex_sig.df)                # 1
table(aldex_sig.df$Comparison)    # Antibiotics vs Control


#### 03. ANCOM-BC2 ####
# Each contrast is a separate two-group model, so normalization and structural
# zero detection are specific to that pair. A single four-group model gives
# different structural zeros and does not reproduce Table 1.
run_ancom_pair <- function(ps.obj, trt, comparison_label){
  keep   <- sample_data(ps.obj)$Treatment %in% c("Control", trt)
  ps_sub <- prune_samples(keep, ps.obj)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

  sd_sub.df <- data.frame(sample_data(ps_sub), check.names = FALSE)
  sd_sub.df$Treatment <- factor(sd_sub.df$Treatment, levels = c("Control", trt))
  sample_data(ps_sub) <- sample_data(sd_sub.df)

  res <- ancombc2(data = ps_sub
                  , fix_formula   = "Treatment"
                  , rand_formula  = NULL
                  , p_adj_method  = "fdr"
                  , prv_cut       = prevalence
                  , lib_cut       = 1000
                  , group         = "Treatment"
                  , struc_zero    = TRUE
                  , neg_lb        = TRUE
                  , alpha         = alpha
                  , iter_control  = list(tol = 1e-5, max_iter = 100, verbose = FALSE)
                  , em_control    = list(tol = 1e-5, max_iter = 100)
                  , lme_control   = lme4::lmerControl()
                  , mdfdr_control = list(fwer_ctrl_method = "holm")
                  , verbose       = FALSE
                  )

  res.df           <- as.data.frame(res$res)
  colnames(res.df) <- make.names(colnames(res.df))
  if("taxon" %in% colnames(res.df)){
    res.df <- rename(res.df, ASV = taxon)
  } else {
    res.df <- rownames_to_column(res.df, var = "ASV")
  }

  trt_std <- make.names(trt)
  res.df %>%
    dplyr::select(ASV
                  , lfc  = !!sym(paste0("lfc_Treatment",  trt_std))
                  , pval = !!sym(paste0("p_Treatment",    trt_std))
                  , padj = !!sym(paste0("q_Treatment",    trt_std))
                  , diff = !!sym(paste0("diff_Treatment", trt_std))
                  ) %>%
    mutate(Comparison = comparison_label)
}

ancom_all.df <- bind_rows(lapply(contrasts.list
                                 , function(x) run_ancom_pair(ps.obj, x$trt, x$label)))
ancom_sig.df <- ancom_all.df %>% filter(diff == TRUE)

nrow(ancom_sig.df)                # 3
table(ancom_sig.df$Comparison)    # 2 Antibiotics, 1 High temperature


#### 04. Table 1 ####
# Manuscript wording for each contrast
label_map <- c("Antibiotics vs Control"      = "Control vs AB"
               , "High temperature vs Control" = "Control vs HT"
               , "Antibiotics + HT vs Control" = "Control vs AB+HT"
               )

# Direction comes from clr_trt_minus_ctrl, never from the sign of diff.btw.
# ALDEx2 sorts conditions alphabetically, so diff.btw is Control-minus-treatment
# for Antibiotics and Antibiotics + HT but treatment-minus-Control for High
# temperature. Keying off it inverts the arrow on two of three contrasts.
#
# The reported value is the ALDEx2 standardised effect size, not a fold change.
# 2^|diff.btw| is a CLR-space quantity and is not comparable to the ANCOM-BC2
# abundance fold changes; putting both in one column invites a false comparison.
table1_aldex.df <- aldex_sig.df %>%
  mutate(Effect        = paste0(ifelse(clr_trt_minus_ctrl > 0, "↑", "↓")
                                , format(round(abs(effect), 2), nsmall = 2))
         , Effect_metric = "ALDEx2 effect size"
         , Comparison    = label_map[Comparison]
         , q_value       = round(we.eBH, 3)
         ) %>%
  dplyr::select(ASV, Comparison, Effect, Effect_metric, q_value) %>%
  mutate(Method = "ALDEx2") %>%
  left_join(tax.df, by = "ASV") %>%
  mutate(across(where(is.character), ~replace_na(.x, "Unclassified")))

# ANCOM-BC2 lfc is Treatment minus Control, so positive means enriched in the
# treatment. Fold change is 2^|lfc|, not exp(|lfc|): base 2 reproduces the
# published values exactly (ASV180 2.6x, ASV35 2.3x, ASV57 4.3x).
table1_ancom.df <- ancom_sig.df %>%
  mutate(fold_change   = round(2^abs(lfc), 1)
         , Effect        = paste0(ifelse(lfc > 0, "↑", "↓"), fold_change, "×")
         , Effect_metric = "fold change (2^LFC)"
         , Comparison    = label_map[Comparison]
         , q_value       = round(padj, 3)
         ) %>%
  dplyr::select(ASV, Comparison, Effect, Effect_metric, q_value) %>%
  mutate(Method = "ANCOM-BC2") %>%
  left_join(tax.df, by = "ASV") %>%
  mutate(across(where(is.character), ~replace_na(.x, "Unclassified")))

table1.df <- bind_rows(table1_aldex.df, table1_ancom.df) %>%
  mutate(Taxonomy = ifelse(Genus != "Unclassified"
                           , paste0(Genus, " (", Family, ")")
                           , paste0(Family, " (unresolved)"))
         ) %>%
  dplyr::select(Comparison, Taxonomy, ASV, Method, Effect, Effect_metric, q_value) %>%
  arrange(Comparison, Method, q_value)

nrow(table1.df)   # 4
table1.df %>% dplyr::select(-ASV)


#### 05. Table S4 ####
# diff.btw and diff.win are kept so the reported effect sizes can be recomputed
# from this file alone. effect is diff.btw / diff.win, so dropping them loses that.
table_s4_aldex.df <- aldex_all.df %>%
  dplyr::select(ASV, Comparison, effect, diff.btw, diff.win, we.ep, we.eBH) %>%
  rename(effect_or_lfc = effect, p_raw = we.ep, p_adj = we.eBH) %>%
  mutate(Method = "ALDEx2") %>%
  left_join(tax.df, by = "ASV")

table_s4_ancom.df <- ancom_all.df %>%
  dplyr::select(ASV, Comparison, lfc, pval, padj) %>%
  rename(effect_or_lfc = lfc, p_raw = pval, p_adj = padj) %>%
  mutate(diff.btw = NA_real_, diff.win = NA_real_, Method = "ANCOMBC2") %>%
  left_join(tax.df, by = "ASV")

table_s4.df <- bind_rows(table_s4_aldex.df, table_s4_ancom.df) %>%
  # Unresolved ranks read as blanks in a spreadsheet unless they are named
  mutate(across(where(is.character), ~replace_na(.x, "Unclassified"))) %>%
  arrange(Method, Comparison, p_adj)


#### 06. Write outputs ####
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

write.csv(table1.df,   file = table1.FN,   row.names = FALSE)
write.csv(table_s4.df, file = table_s4.FN, row.names = FALSE)

# Go to scripts/05_indicator_species.R
