## =======================================================================
## Script Name: run_metabolism5_GWAS.R
## Description: Run multivariate GWAS for the shared metabolic liability 
##              (Single Common Factor Model) using Genomic SEM on an HPC.
## =======================================================================

# 1. Set the working directory (Change this to your actual server path)
setwd("mnt/metabolism5_GWAS")

# 2. Load necessary R packages
require(GenomicSEM)
require(data.table)

starttime <- Sys.time()
print(paste("Start time:", starttime))

# 3. Load the LDSC output covariance matrix
print("Loading LDSC covariance matrix...")
load("metabolism5_S_LDSCOutput.RData") 

print("Loading prepped summary statistics...")
print(Sys.time() - starttime)

# 4. Load the prepped summary statistics containing all genome-wide SNPs
sumstats <- fread("prepped_sumstats_metabolism5.txt")

print("Data loaded successfully.")
print(Sys.time() - starttime)

# 5. Get environment variables from SLURM array task ID
args <- commandArgs(trailingOnly = TRUE)
job_id <- as.numeric(args[1])

# Fallback for local testing if job_id is missing
if (is.na(job_id)) {
  job_id <- 1
}

# 6. Split millions of SNPs into manageable chunks for parallel processing
chunk_size <- 10000 
snp_index <- 1:nrow(sumstats)
chunks <- split(snp_index, ceiling(seq_along(snp_index) / chunk_size))
chunk_idx <- job_id
current_chunk <- chunks[[chunk_idx]]

# 7. Define the single common factor model for the 5 metabolic traits
# F1 represents the shared metabolic liability
model <- "
F1 =~ NA*GOUT + HYPTENSESS + DISLIPIDIMIA + NAFLD + T2D
F1 ~~ 1*F1
GOUT ~~ a*GOUT
a > .0001
HYPTENSESS ~~ b*HYPTENSESS
b > .0001
DISLIPIDIMIA ~~ c*DISLIPIDIMIA
c > .0001
NAFLD ~~ d*NAFLD
d > .0001
T2D ~~ e*T2D
e > .0001
F1 ~ SNP
"

# 8. Specify the parameter to extract (SNP effect on the common factor F1)
sub_extract <- c("F1~SNP")

print(paste("Running multivariate GWAS for chunk:", chunk_idx))
print(Sys.time() - starttime)

# 9. Run the Genomic SEM userGWAS function
results <- userGWAS(
  covstruc = LDSCoutput, 
  SNPs = sumstats[current_chunk, ], 
  model = model, 
  parallel = TRUE, 
  sub = sub_extract, 
  cores = 36,                 
  GC = "none", 
  smooth_check = TRUE, 
  printwarn = TRUE,
  toler = 1e-30,              
  fix_measurement = TRUE,     
  Q_SNP = TRUE                
)

print("Writing GWAS output to file...")
print(Sys.time() - starttime)

# 10. Save the chunk results to a tab-separated file
output_filename <- paste0("metabolism5_factor_F1_chunk_", chunk_idx, ".tab")

write.table(results[[1]], 
            file = output_filename, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

print(paste("Successfully finished chunk:", chunk_idx))
print(Sys.time() - starttime)
