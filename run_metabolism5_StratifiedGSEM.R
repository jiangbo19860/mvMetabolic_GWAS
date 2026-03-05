#!/usr/bin/env Rscript

# ====================================================================================
# Script Name: run_metabolism5_StratifiedGSEM.R
# Description: Stratified Genomic Structural Equation Modeling (S-GSEM) for 
#              the shared genetic architecture of 5 metabolic conditions.
# ====================================================================================

# 1. Load required packages
require(stringr)
require(GenomicSEM)

# 2. Set working directory (Change this to your actual server path)
setwd("mnt/metabolism5_SGSEM/")

print("Starting Stratified Genomic SEM for 5 Metabolic Diseases...")

# 3. Define the input summary statistics files
# Ensure these are the correctly munged sumstats files located in your working directory
sumstats_files <- c(
  "GOUT.sumstats.gz",
  "HYPTENSESS.sumstats.gz",
  "LIPOPROT.sumstats.gz",
  "NAFLD.sumstats.gz",
  "T2D.sumstats.gz"
)

# 4. Define trait names
trait_names <- c("GOUT", "HYPTENSESS", "LIPOPROT", "NAFLD", "T2D")

# 5. Define sample prevalences (proportion of cases in the GWAS samples)
# Assuming 0.5 for all as effective sample size (Neff) was used for binary traits
sample_prev <- c(0.5, 0.5, 0.5, 0.5, 0.5)

# 6. Define population prevalences (based on large-scale epidemiological data)
population_prev <- c(0.0066, 0.33, 0.241, 0.3, 0.1111)

# 7. Define the Single Common Factor Model for Metabolic Liability
# F1 represents the shared metabolic liability. 
# Residual variances are constrained to be > 0.0001 to prevent negative variances.
model <- "
F1 =~ NA*GOUT + HYPTENSESS + LIPOPROT + NAFLD + T2D
F1 ~~ 1*F1
GOUT ~~ a*GOUT
a > .0001
HYPTENSESS ~~ b*HYPTENSESS
b > .0001
LIPOPROT ~~ c*LIPOPROT
c > .0001
NAFLD ~~ d*NAFLD
d > .0001
T2D ~~ e*T2D
e > .0001
"

print("Model specified successfully. Initiating Stratified Genomic SEM...")

# 8. Run Stratified Genomic SEM
# Note: Ensure that the 'ld', 'wld', and 'frq' paths correctly point to your local 
# downloaded 1000 Genomes Phase 3 baseline reference folders.
S_gSEM_results <- S_gSEM(
  Sumstatsfile = sumstats_files,
  trait.names = trait_names,
  sample.prev = sample_prev,
  population.prev = population_prev,
  ld = "LDSCORE_1000G_Phase3_baselineLD_v2.2_ldscores",  # Path to baseline LD scores
  wld = "weights_hm3_no_hla",                             # Path to LD weights
  frq = "1000G_Phase3_frq",                               # Path to allele frequencies
  model = model,
  params = c("F1~~F1"),                                   # Extract heritability of the latent factor
  fix = "regressions",
  std.lv = FALSE,
  rm_flank = TRUE,
  tau = FALSE,
  base = TRUE,
  toler = NULL,
  fixparam = NULL,
  save_name = "metabolism5_SGSEM_enrichment",
  save_path = "./"
)

print("Stratified Genomic SEM completed successfully!")
print("Results are saved in the current directory as 'metabolism5_SGSEM_enrichment.RData' and CSV files.")
