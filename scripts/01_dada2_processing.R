# Build the phyloseq object the rest of the pipeline reads, from raw amplicon reads
#  requires paired-end fastq in data/raw/ and the SILVA v138 training set
#  initialized 2026-09-07
#  James A. Dennis-Orr (james.arnold.do@gmail.com)
#
#  Raw reads are not in this repository. They go to NCBI BioProject PRJXXXXXX on
#  acceptance. The script is here so the provenance of the deposited object is
#  visible; to reproduce the manuscript, start from scripts/02.
#
#  It writes the object unfiltered apart from DADA2 quality control: algae
#  samples, zero-count taxa, prevalence and rarefaction are all handled by the
#  downstream scripts, so every ecological choice stays where it is applied.
#
#  Adapted from the DADA2 big-data workflow (Callahan et al. 2016,
#  Nat Methods 13:581-583). Taxonomy from SILVA v138 (Quast et al. 2013).

#### 00. Front Matter ####
# Clear space
# rm(list=ls())

# Load libraries
#BiocManager::install("dada2")
library(dada2)
#BiocManager::install("phyloseq")
library(phyloseq)
#BiocManager::install("DECIPHER")
library(DECIPHER)
#install.packages("phangorn")
library(phangorn)

# User set variables
raw_dir.FN   <- "data/raw"
out_dir.FN   <- "data/processed"
silva.FN     <- "data/raw/silva_nr99_v138_train_set.fa.gz"   # zenodo.org/record/3986799
metadata.FN  <- "data/processed/sample_metadata.csv"
out.FN       <- "data/processed/James_MU42022_filtered_not_normalized.RDS"
asv_seqs.FN  <- "data/processed/ASV_sequences.txt"
trunc_len    <- c(260, 200)   # set these from plotQualityProfile(), see below
max_ee       <- c(2, 5)
seed         <- 100

set.seed(seed)


#### 01. Locate reads ####
fnFs <- sort(list.files(raw_dir.FN, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(raw_dir.FN, pattern = "_R2_001.fastq", full.names = TRUE))
sample_names <- sapply(strsplit(basename(fnFs), split = "_"), `[`, 1)

length(fnFs)   # 67 libraries
stopifnot(length(fnFs) == length(fnRs))


#### 02. Filter and trim ####
# Inspect quality before committing to trunc_len; the defaults above suit this run
# plotQualityProfile(fnFs[1:4])
# plotQualityProfile(fnRs[1:4])

filt_dir.FN <- file.path(raw_dir.FN, "filtered")
dir.create(filt_dir.FN, showWarnings = FALSE, recursive = TRUE)

filtFs <- file.path(filt_dir.FN, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir.FN, paste0(sample_names, "_R_filt.fastq.gz"))

filt_out <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs
                          , truncLen = trunc_len
                          , maxEE    = max_ee
                          , truncQ   = 2
                          , rm.phix  = TRUE
                          , compress = TRUE
                          , multithread = TRUE
                          )
head(filt_out)


#### 03. Infer ASVs ####
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

derepFs <- derepFastq(filtFs); names(derepFs) <- sample_names
derepRs <- derepFastq(filtRs); names(derepRs) <- sample_names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

mergers   <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtab    <- makeSequenceTable(mergers)
seqtab.nc <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

round(100 * sum(seqtab.nc) / sum(seqtab), 1)   # percent of reads kept through chimera removal
table(nchar(getSequences(seqtab.nc)))          # amplicon length distribution


#### 04. Assign taxonomy ####
tax.mat <- assignTaxonomy(seqtab.nc, refFasta = silva.FN, multithread = TRUE)

# Short names in place of the sequences themselves, with the mapping kept
asv_seqs  <- colnames(seqtab.nc)
asv_names <- paste0("ASV", seq_along(asv_seqs))
colnames(seqtab.nc) <- asv_names
rownames(tax.mat)   <- asv_names

write.table(data.frame(ASV = asv_names, Sequence = asv_seqs)
            , file = asv_seqs.FN, sep = "\t", quote = FALSE, row.names = FALSE)


#### 05. Phylogenetic tree ####
# UniFrac needs a tree. Note it comes out unrooted; script 03 midpoint-roots it
# rather than letting phyloseq pick a root at random on every call.
seqs <- getSequences(seqtab.nc); names(seqs) <- asv_names
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = FALSE)

phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
tree_nj     <- NJ(dist.ml(phang.align))
fit_gtr     <- optim.pml(update(pml(tree_nj, data = phang.align), k = 4, inv = 0.2)
                         , model         = "GTR"
                         , optInv        = TRUE
                         , optGamma      = TRUE
                         , rearrangement = "stochastic"
                         , control       = pml.control(trace = 0)
                         )


#### 06. Assemble and write ####
metadata.df <- read.csv(file = metadata.FN, row.names = 1)

ps.obj <- phyloseq(otu_table(seqtab.nc, taxa_are_rows = FALSE)
                   , tax_table(tax.mat)
                   , sample_data(metadata.df)
                   , phy_tree(fit_gtr$tree)
                   )

ntaxa(ps.obj)             # 7326
nsamples(ps.obj)          # 64, algae included; scripts 02 onward drop them
min(sample_sums(ps.obj))  # 1682, not rarefied

dir.create(out_dir.FN, showWarnings = FALSE, recursive = TRUE)
saveRDS(ps.obj, file = out.FN)

# Go to scripts/02_alpha_diversity.R
