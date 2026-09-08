# Model Vibrio challenge survival by early-life treatment, and draw Figure 5
#  requires data/challenge/vibrio_survival.csv
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)
#
#  Challenge well enters as a random effect. Tank is perfectly nested in
#  treatment (one rearing tank each), so no term can separate the two.

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
#install.packages("here")
library(here)
#install.packages("survival")
library(survival)
#install.packages("coxme")
library(coxme)
#install.packages("tidyverse")
library(tidyverse)

source(here::here("scripts/00_colors.R"))

# User set variables
input.FN     <- "data/challenge/vibrio_survival.csv"
cox.FN       <- "tables/cox_results.csv"
ph_test.FN   <- "tables/cox_ph_test.csv"
model.FN     <- "tables/coxme_model_summary.txt"
fig5.FN      <- "figures/Figure5_cox_forest"
per_well     <- 10       # spat per challenge well
seed         <- 100

set.seed(seed)


#### 01. Load ####
dat.df <- read.csv(file = input.FN)
stopifnot(all(c("TE", "Binary", "treatment", "animal") %in% names(dat.df)))

nrow(dat.df)              # 160
table(dat.df$treatment)   # 40 each
sort(unique(dat.df$TE))   # 3 4 5, only two of which carry deaths


#### 02. Reconstruct well assignments ####
# Animals were entered sequentially by well, ten per well, four wells per
# treatment. Confirmed against the lab book. Positions 1-10 are well 1, and so on.
dat.df <- dat.df %>%
  arrange(animal) %>%
  group_by(treatment) %>%
  mutate(rel_pos  = row_number()
         , well_num = ceiling(rel_pos / per_well)
         , well     = paste0(treatment, "_W", well_num)
         ) %>%
  ungroup()

length(unique(dat.df$well))   # 16

dat.df <- dat.df %>%
  mutate(treatment = factor(treatment
                            , levels = c("A", "B", "C", "D")
                            , labels = TREATMENT_ORDER)
         , well      = factor(well)
         )

# Mortality quoted in Results 3.5
dat.df %>%
  group_by(treatment) %>%
  summarise(n = n(), deaths = sum(Binary), pct = round(100 * mean(Binary), 1))
# Control 15.0, Antibiotics 47.5, High temperature 70.0, Antibiotics + HT 67.5


#### 03. Mixed-effects Cox model ####
coxme.fit <- coxme(Surv(TE, Binary) ~ treatment + (1 | well), data = dat.df)
coxme.fit
VarCorr(coxme.fit)


#### 04. Proportional hazards check ####
# cox.zph has no coxme method. If PH holds for the fixed effects alone it holds
# for the same fixed effects with a well-level random intercept added.
cox_fixed.fit <- coxph(Surv(TE, Binary) ~ treatment, data = dat.df, ties = "efron")
zph           <- cox.zph(cox_fixed.fit)
zph   # global chisq 3.18, df 3, p 0.36

ph_test.df <- as.data.frame(zph$table) %>%
  rownames_to_column(var = "term") %>%
  rename(p_value = p)


#### 05. Hazard ratios ####
# coxme has no broom::tidy method, so the table is built from the coefficient
# vector and the diagonal of the variance-covariance matrix.
fixed_coef <- fixef(coxme.fit)
fixed_se   <- sqrt(diag(vcov(coxme.fit)))

hr.df <- data.frame(Treatment = names(fixed_coef)
                    , log_HR    = fixed_coef
                    , SE        = fixed_se
                    ) %>%
  mutate(HR        = exp(log_HR)
         , CI_lower  = exp(log_HR - 1.96 * SE)
         , CI_upper  = exp(log_HR + 1.96 * SE)
         , p_value   = 2 * pnorm(-abs(log_HR / SE))
         , Treatment = gsub("^treatment", "", Treatment)
         , Treatment = factor(Treatment
                              , levels = rev(c("Antibiotics", "High temperature", "Antibiotics + HT")))
         )

hr.df %>% dplyr::select(Treatment, HR, CI_lower, CI_upper, p_value)
# HT 7.96 (1.94-32.65) p 0.004 | AB+HT 7.40 (1.79-30.48) p 0.006 | AB 3.66 (0.86-15.60) p 0.079

cox.df <- hr.df %>%
  dplyr::select(Treatment, HR, CI_lower, CI_upper, p_value) %>%
  mutate(across(c(HR, CI_lower, CI_upper), ~round(.x, 2))
         , p_value = round(p_value, 4)
         )


#### 06. Figure 5 ####
hr.df <- hr.df %>%
  mutate(label = sprintf("%.2f (%.2f–%.2f)\np = %.3f", HR, CI_lower, CI_upper, p_value))

fig5.plot <- ggplot(hr.df, aes(x = HR, y = Treatment, xmin = CI_lower, xmax = CI_upper)) +
  geom_errorbarh(height = 0.18, linewidth = 0.8, colour = "grey40") +
  geom_point(size = 3.5, aes(colour = Treatment)) +
  geom_text(aes(label = label), nudge_y = -0.30, hjust = 0.5
            , size = 3.2, colour = "black", lineheight = 0.9) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6, colour = "grey30") +
  # Log scale because the intervals span nearly two orders of magnitude
  scale_x_log10(name = "Hazard ratio (log scale)"
                , breaks = c(1, 2, 5, 10, 20)
                , limits = c(0.5, 30)
                ) +
  scale_colour_manual(values = TREATMENT_COLORS, guide = "none") +
  theme_classic(base_size = 14) +
  theme(axis.title.y = element_blank()
        , axis.ticks.y = element_blank()
        , axis.line.y  = element_blank()
        , axis.text.y  = element_text(size = 13)
        , plot.title   = element_text(face = "bold", size = 14, hjust = 0.5)
        , plot.margin  = margin(10, 20, 30, 10)
        ) +
  labs(title = "Hazard Ratios Relative to Control\n(coxme, random effect = challenge well, 95% CI)")


#### 07. Write outputs ####
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables",  showWarnings = FALSE, recursive = TRUE)

write.csv(cox.df,     file = cox.FN,     row.names = FALSE)
write.csv(ph_test.df, file = ph_test.FN, row.names = FALSE)

# Each section is labelled so the saved file can be read on its own
model_sections.list <- list("coxme model"                      = coxme.fit
                            , "Random effect variance"           = VarCorr(coxme.fit)
                            , "Standard Cox (PH diagnostics only)" = summary(cox_fixed.fit)
                            , "PH assumption (Schoenfeld)"       = zph
                            , "HR table"                         = hr.df %>%
                                dplyr::select(Treatment, HR, CI_lower, CI_upper, p_value)
                            )

sink(file = model.FN)
for(nm in names(model_sections.list)){
  writeLines(c("", paste("===", nm, "==="), ""))
  print(model_sections.list[[nm]])
}
sink()

ggsave(filename = paste0(fig5.FN, ".png"), plot = fig5.plot, width = 6, height = 4, dpi = 300)
ggsave(filename = paste0(fig5.FN, ".pdf"), plot = fig5.plot, width = 6, height = 4, device = cairo_pdf)

# Go to scripts/08_dataset_summary.R
