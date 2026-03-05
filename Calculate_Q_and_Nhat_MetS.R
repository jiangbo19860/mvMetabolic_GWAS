#!/usr/bin/env Rscript

# ====================================================================================
# Script Name: Calculate_Q_and_Nhat_MetS.R
# Description: Post-GWAS Quality Control, Heterogeneity (Q_SNP) assessment, 
#              and Effective Sample Size (N_hat) calculation for the shared 
#              metabolic liability factor (F1).
# ====================================================================================

# 1. Load required packages
require(data.table)
require(dplyr)
require(R.utils)

# Set working directory (Change to your actual path)
setwd("mnt/metabolism5_GWAS/")

print("Loading multivariate GWAS summary statistics...")

# 2. Load the concatenated GWAS results generated from the previous step
Data <- fread("Final_Multivariate_GWAS_Metabolism.txt", header = TRUE)

initial_snp_count <- nrow(Data)
print(paste("Initial number of SNPs:", initial_snp_count))

# ====================================================================================
# Step 1: Quality Control (Remove warnings and smoothed SNPs)
# ====================================================================================

# Remove SNPs that produced warnings during the Genomic SEM estimation
Data <- subset(Data, warning == 0)
print(paste("SNPs remaining after removing warnings:", nrow(Data)))

# Remove SNPs that required smoothing (Z_smooth != 0) to ensure high-confidence estimates
# This strictly retains SNPs where the genetic covariance matrix was naturally positive definite
Data <- subset(Data, Z_smooth == 0)
print(paste("Number of final high-quality SNPs:", nrow(Data)))

# ====================================================================================
# Step 2: Assess Heterogeneity (Q_SNP)
# ====================================================================================

# Define Q_SNP significance threshold for a single common factor model
# (Standard genome-wide significance threshold)
Qsig <- 5e-8

# Count how many SNPs show significant heterogeneity 
# (i.e., their effects are NOT fully mediated by the shared metabolic factor)
num_Qsig <- sum(Data$Q_pval <= Qsig, na.rm = TRUE)
print(paste("Number of SNPs with significant heterogeneity (Q_pval <= 5e-8):", num_Qsig))

# Count how many genome-wide significant SNPs are ALSO highly heterogeneous
sig_SNPs <- subset(Data, Pval_Estimate <= 5e-8)
num_sig_and_Qsig <- sum(sig_SNPs$Q_pval <= Qsig, na.rm = TRUE)
print(paste("Number of GW-significant SNPs that are also heterogeneous:", num_sig_and_Qsig))

# Save the full summary statistics (including Q-significant SNPs) for archival
print("Saving full summary statistics (with Q SNPs)...")
write.table(Data, file = "metabolism5_factor_F1_sumstats_full.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
gzip("metabolism5_factor_F1_sumstats_full.txt", overwrite = TRUE)

# Extract and save positional information of top Q-significant SNPs for potential future exploration
Q_SNPs_list <- subset(Data, Q_pval <= Qsig, select = c(SNP, CHR, BP))
write.table(Q_SNPs_list, file = "metabolism5_factor_F1_topQsnps.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")

# ====================================================================================
# Step 3: Remove Q-significant SNPs for downstream analyses
# ====================================================================================

# For downstream analyses (like MR and MTWAS) that require the genetic instrument 
# to act strictly through the shared metabolic liability, we exclude heterogeneous SNPs.
Data_noQ <- subset(Data, Q_pval > Qsig | is.na(Q_pval))

print(paste("SNPs remaining after excluding Q-significant SNPs:", nrow(Data_noQ)))

# ====================================================================================
# Step 4: Calculate Effective Sample Size (N_hat)
# ====================================================================================

# N_hat is calculated using well-imputed, common variants (MAF between 10% and 40%)
# Formula: N_hat = mean( 1 / (2 * MAF * (1 - MAF) * SE^2) )
CorrelatedFactors <- subset(Data_noQ, MAF <= 0.4 & MAF >= 0.1)

# Ensure the standard error column is correctly referenced ('se_c' from Genomic SEM)
N_hat <- mean(1 / ((2 * CorrelatedFactors$MAF * (1 - CorrelatedFactors$MAF)) * (CorrelatedFactors$se_c^2)), na.rm = TRUE)

print(paste("Calculated Effective Sample Size (N_hat) for the shared factor:", round(N_hat)))

# Assign the calculated N_hat to all SNPs in the dataframe
Data_noQ$N <- round(N_hat)

# Clean up memory
rm(CorrelatedFactors)

# ====================================================================================
# Step 5: Final Output
# ====================================================================================

# Assess the number of genome-wide significant SNPs in the strict, homogeneous dataset
final_sig_count <- sum(Data_noQ$Pval_Estimate <= 5e-8, na.rm = TRUE)
print(paste("Final number of homogeneous, genome-wide significant SNPs:", final_sig_count))

print("Saving final strict summary statistics (without Q SNPs)...")

# Write the final summary statistics file to be used for FUMA, MTWAS, and MR
write.table(Data_noQ, file = "metabolism5_factor_F1_sumstats_noQ.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")
gzip("metabolism5_factor_F1_sumstats_noQ.txt", overwrite = TRUE)

print("Post-GWAS processing completed successfully!")
