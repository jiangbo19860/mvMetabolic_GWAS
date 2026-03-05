## Run GWAS of the single common factor model on the Computing Cluster

setwd("/mnt/metabolism5_GWAS")

# Load required R packages
require(gdata)
require(lavaan)
require(doParallel)
require(stringr)
require(Matrix)
require(R.utils)
require(dplyr)
require(utils)
require(GenomicSEM)
require(data.table)

starttime <- Sys.time()
cat("Start time:", format(starttime), "\n")

cat("Loading data...\n")

# Load LDSC results object
# The RData should contain an object named LDSCoutput
load("metabolism5_S_gSEM_S_LDSCOutput.RData") 

cat("Processing data...\n")

# Load and read summary statistics file containing SNP information
# (should be a file processed and merged by munge function with genome-wide SNPs)
sumstats <- fread("prepped_sumstats_metabolism5.txt")

cat("Data loaded. Elapsed time:", format(Sys.time() - starttime), "\n")

# Get environment variable (Task ID) from SLURM array task for chunk specification
args <- commandArgs(trailingOnly = TRUE)
job_id <- if (length(args) > 0) as.numeric(args[1]) else 1

# Divide millions of SNPs into chunks (adjust chunk_size based on node configuration)
chunk_size <- 10000 
snpj <- 1:nrow(sumstats)
chunks <- split(snpj, ceiling(seq_along(snpj) / chunk_size))
chunk_idx <- min(job_id, length(chunks))
chunk <- chunks[[chunk_idx]]

# Define single factor structural equation model for metabolism 5 diseases
# Note: SNP effect estimation requires explicit "F1 ~ SNP" line for userGWAS
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
F1 ~ SNP
"

# Specify regression paths to extract (extract only F1 explained by SNP)
sub <- c("F1~SNP")

cat("Running GWAS for chunk:", chunk_idx, "\n")
cat("Elapsed time:", format(Sys.time() - starttime), "\n")

# Run multivariate GWAS
# toler=1e-30 handles singular matrix errors; fix_measurement=TRUE speeds up computation
results <- userGWAS(
    covstruc = LDSCoutput, 
    SNPs = sumstats[chunk, ], 
    model = model, 
    parallel = TRUE, 
    sub = sub, 
    cores = 36,
    GC = "none", 
    smooth_check = TRUE, 
    printwarn = TRUE,
    toler = 1e-30,
    fix_measurement = TRUE,
    Q_SNP = TRUE
)

cat("Writing GWAS output...\n")

# Save output results
write.table(results[[1]], 
                        file = paste0("metabolism5_factor_F1_", chunk_idx, ".tab"), 
                        sep = "\t", 
                        row.names = FALSE, 
                        quote = FALSE)

cat("Finished chunk:", chunk_idx, "\n")
cat("Total elapsed time:", format(Sys.time() - starttime), "\n")
