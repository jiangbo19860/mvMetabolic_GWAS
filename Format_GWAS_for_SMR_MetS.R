#!/usr/bin/env Rscript

# ====================================================================================
# Script Name: Format_GWAS_for_SMR_MetS.R
# Description: Align the metabolic liability GWAS SNPs with the A1 and A2 alleles 
#              of the 1000G reference panel. Convert MAF to EAF (Effect Allele 
#              Frequency) using major allele information from dbSNP, and format 
#              the output strictly for SMR (.ma format).
# ====================================================================================

# 1. Load required packages
library(dplyr)
library(data.table)
library(R.utils)

# Set working directory (Change to your actual SMR directory path)
setwd("/your/server/path/metabolism5_SMR/")

print("Starting GWAS summary statistics formatting for SMR...")

# ====================================================================================
# Step 1: Prepare the 1000G reference panel combined with dbSNP major alleles
# ====================================================================================
# Note: This assumes you have already generated the '1000G_EUR_Phase3_bim_all_chr_with_dbSNP_major_allele.txt.gz'
# as per the standard SMR preparation protocol.

print("Loading 1000G reference panel with dbSNP major allele information...")
combined_1000G_major <- fread("1000G_EUR_Phase3_bim_all_chr_with_dbSNP_major_allele.txt.gz")
print(paste("Reference panel loaded. Number of SNPs:", nrow(combined_1000G_major)))


# ====================================================================================
# Step 2: Load and align the Metabolic Liability (F1) GWAS summary statistics
# ====================================================================================

print("Loading post-QC, Q-filtered metabolic F1 GWAS summary statistics...")
# Read the file generated from the previous QC step
F1 <- fread("../metabolism5_GWAS/metabolism5_factor_F1_sumstats_noQ.txt.gz")

print(paste("Initial number of F1 SNPs:", nrow(F1)))

# Check for duplicated SNPs in the GWAS file
duplicate_count <- sum(duplicated(F1$SNP))
print(paste("Number of duplicated SNPs in F1:", duplicate_count))
if(duplicate_count > 0) {
  F1 <- subset(F1, !duplicated(F1$SNP))
}

# Merge the F1 GWAS with the 1000G reference panel by rsID
print("Merging F1 GWAS with 1000G reference panel...")
F1 <- merge(F1, combined_1000G_major, by = "SNP")
print(paste("SNPs remaining after merge:", nrow(F1)))


# ====================================================================================
# Step 3: Allele Flipping and Mismatch Removal
# ====================================================================================

# Flip the effect size (est) if the effect allele in GWAS (A1.x) matches 
# the non-effect allele in the reference (A2.y)
print("Flipping beta estimates where alleles are reversed...")
F1$est <- ifelse(F1$A1.x != F1$A1.y & F1$A1.x == F1$A2.y, F1$est * -1, F1$est)

# Remove SNPs that do not match either A1 or A2 in the reference file
print("Removing SNPs with unresolvable allele mismatches...")
F1 <- subset(F1, !(F1$A1.x != F1$A1.y & F1$A1.x != F1$A2.y))
F1 <- subset(F1, !(F1$A2.x != F1$A2.y & F1$A2.x != F1$A1.y))
print(paste("SNPs remaining after mismatch removal:", nrow(F1)))


# ====================================================================================
# Step 4: Calculate Effect Allele Frequency (EAF)
# ====================================================================================
# SMR strictly requires EAF. We calculate EAF using the Major Allele from dbSNP.
# We use A1.y (from the reference) as the new effect allele since betas were flipped to match it.

print("Calculating Effect Allele Frequency (EAF)...")
F1$EAF <- ifelse(F1$A1.y != F1$MajorAllele & F1$A2.y == F1$MajorAllele, F1$MAF, 1 - F1$MAF)

# Check EAF ranges to ensure calculations are valid
print("EAF calculation complete. EAF range:")
print(range(F1$EAF, na.rm = TRUE))


# ====================================================================================
# Step 5: Format and Save the Final .ma File for SMR
# ====================================================================================
# The SMR software strictly requires exactly these 8 columns in this specific order:
# SNP, A1, A2, freq, b, se, p, n

print("Formatting columns to SMR standard (.ma)...")
F1_SMR <- F1 %>% select(SNP, A1.y, A2.y, EAF, est, se_c, Pval_Estimate, N)

# Rename columns to match SMR requirements
names(F1_SMR) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")

output_file <- "metabolism5_factor_F1_SMR.ma"
print(paste("Writing final SMR file to:", output_file))

write.table(F1_SMR, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

print("SMR formatting successfully completed!")
