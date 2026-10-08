# Early-Life Microbiome Perturbations Shape Disease Susceptibility in Pacific Oyster (*Crassostrea gigas*)

**Authors:** James A. Dennis-Orr\*, Marissa D. Wright-LaGreca, Timothy J. Green, Andrew H. Loudon

**Target journal:** Applied and Environmental Microbiology

A 24-hour antibiotic or thermal perturbation applied on Day 1 post-fertilization changes which
bacteria colonize Pacific oyster larvae. Those animals die at higher rates when challenged with
*Vibrio aestuarianus* three months later. Alpha diversity never separates the treatments; community
composition does, and 90 ASVs remain treatment-specific indicators at Day 15.

---

## Data

**Raw sequence reads:** NCBI SRA, BioProject PRJXXXXXX, deposited upon acceptance.

**Processed data:** included here. `data/` holds everything scripts 02 through 08 need.

```
data/
├── processed/
│   ├── James_MU42022_filtered_not_normalized.RDS   phyloseq object, unrarefied counts
│   └── sample_metadata.csv                         sample metadata
└── challenge/
    └── vibrio_survival.csv                         challenge survival, 160 animals (2 low-dose controls excluded in 07)
```

The RDS carries the OTU table, taxonomy (SILVA v138), sample data and the phylogenetic tree.
Scripts apply their own filters to it rather than reading a pre-filtered object, so every
exclusion is visible in the code.

---

## Running the analyses

Open the scripts in RStudio and step through them, or source them all:

```r
source("run_all.R")
```

Run from the repository root. All paths are relative to it; there is no `setwd()` anywhere.
Script 01 is skipped by `run_all.R` because it needs the raw reads.

| Script | Produces |
|--------|----------|
| `00_colors.R` | shared colour, shape and linetype scales, sourced by 02 to 07 |
| `01_dada2_processing.R` | the RDS in `data/processed/`. Needs raw reads, **not runnable here** |
| `02_alpha_diversity.R` | Figure 1, Table S1 |
| `03_beta_diversity.R` | Figure 2, Figure S1, Table S3 |
| `04_differential_abundance.R` | Table 1, Table S4 |
| `05_indicator_species.R` | Figure 4, Table 2, Table S2, two IndVal detail tables |
| `06_composition_barplots.R` | Figure 3 |
| `07_survival_analysis.R` | Figure 5, Cox tables, full model printout |
| `08_dataset_summary.R` | prints the dataset numbers quoted in Methods 2.5 and 2.6 |

Figures are written to `figures/` as PNG at 300 dpi and PDF. Tables go to `tables/`.

Figure 6 in the manuscript comes from a separate follow-up challenge conducted in 2025. Its data
and analysis are not part of this repository.

Script 04 is the slow one. ALDEx2 draws 128 Monte Carlo instances for each of three comparisons
and ANCOM-BC2 fits three separate models, so allow several minutes.

---

## Packages

R 4.4 or later.

```r
install.packages(c("here", "tidyverse", "vegan", "rstatix", "patchwork",
                   "ggrepel", "survival", "coxme", "indicspecies", "phangorn"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "metagenomeSeq", "ALDEx2", "ANCOMBC"))
```

Script 01 additionally needs `dada2` and `DECIPHER`, plus the SILVA v138 training set
(`silva_nr99_v138_train_set.fa.gz`, from <https://zenodo.org/record/3986799>) placed in `data/raw/`
alongside the fastq files.

Each script carries its install command commented directly above the matching `library()` call.

---

## Parameters

Everything settable sits in a `# User set variables` block at the top of each script.

| Parameter | Value |
|-----------|-------|
| Rarefaction depth (alpha diversity) | 5,444 reads |
| Beta diversity normalization | CSS (metagenomeSeq) |
| Prevalence filter | present in ≥10% of samples, giving 739 ASVs |
| Matched sampling days | Days 1, 3, 6, 15 |
| n, alpha diversity | 44 (one sample falls below the rarefaction depth) |
| n, beta diversity | 45 (CSS does not discard samples) |
| PERMANOVA and IndVal permutations | 9,999 |
| ALDEx2 Monte Carlo instances | 128 |
| Multiple testing, differential abundance | BH (FDR), with Holm across directional tests |
| Multiple testing, pairwise PERMANOVA and Dunn | Bonferroni |
| Random seed | 100 |

Two details worth knowing before modifying the code:

**The phylogenetic tree is unrooted.** UniFrac needs a root, and phyloseq assigns one at random
when none is given, which makes the distances differ on every call. Script 03 midpoint-roots the
tree first. Removing that line makes the UniFrac results irreproducible.

**`MASS::select` masks `dplyr::select`.** `metagenomeSeq`, `ALDEx2` and `ANCOMBC` all pull in MASS,
and attached packages accumulate across a `run_all.R` session. Every script therefore writes
`dplyr::select()` in full. Shortening it breaks the run.

The unweighted UniFrac PCoA is deliberately not drawn. Its dispersion is unequal between treatments
(betadisper F = 4.74, p = 0.0056), so a PERMANOVA on it would not be interpretable. It appears in
Table S3 only.

---

## License

MIT, see `LICENSE`.
