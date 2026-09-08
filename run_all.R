# Source every analysis script in order and report what each produced
#  requires the processed data in data/, see README
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)
#
#  Script 01 is skipped: it needs the raw reads, which are on the SRA rather
#  than in this repository.
#
#  Each script is sourced into its own environment. That keeps their objects
#  apart, but attached packages still accumulate across the run, so every
#  script qualifies dplyr::select() rather than relying on the search path.

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# User set variables
scripts.FN <- c("scripts/02_alpha_diversity.R"
                , "scripts/03_beta_diversity.R"
                , "scripts/04_differential_abundance.R"
                , "scripts/05_indicator_species.R"
                , "scripts/06_composition_barplots.R"
                , "scripts/07_survival_analysis.R"
                , "scripts/08_dataset_summary.R"
                )


#### 01. Run ####
results.df <- data.frame(script  = basename(scripts.FN)
                         , status  = NA_character_
                         , elapsed = NA_real_
                         , stringsAsFactors = FALSE
                         )

for(i in seq_along(scripts.FN)){
  message("[", i, "/", length(scripts.FN), "] ", basename(scripts.FN[i]))
  t0 <- proc.time()["elapsed"]

  tryCatch({
    source(scripts.FN[i], local = new.env())
    results.df$status[i]  <- "OK"
    results.df$elapsed[i] <- round(proc.time()["elapsed"] - t0, 1)
  }, error = function(e){
    results.df$status[i]  <<- paste("ERROR:", conditionMessage(e))
    results.df$elapsed[i] <<- round(proc.time()["elapsed"] - t0, 1)
  })
}


#### 02. Report ####
results.df

list.files("figures", pattern = "\\.(png|pdf)$")   # 12 files, six figures in two formats
list.files("tables",  pattern = "\\.(csv|txt)$")   # 11 files

if(any(results.df$status != "OK")){
  warning("Not every script completed, see the status column")
}

# End
