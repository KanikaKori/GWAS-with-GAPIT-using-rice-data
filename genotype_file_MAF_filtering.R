# 8. Filtering for minor allele frequency (MAF ≥ 0.05) and missing data rate (≤ 10%)
# Load necessary library
library(dplyr)

# Step A: Load your hapmap file 
hapmap <- read.delim("genotype_file.hmp.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Step B: Extract only genotype columns (exclude first 11 metadata columns) 
geno_data <- hapmap[, -(1:11)]

# Step C: Define a function to calculate MAF
calculate_maf <- function(genotypes) {
  # Remove missing data
  alleles <- unlist(strsplit(na.omit(genotypes), split = ""))
  alleles <- alleles[alleles %in% c("A", "T", "C", "G")]
  
  # If no valid alleles, return NA
  if (length(alleles) < 2) return(NA)
  
  # Count allele frequencies
  allele_freqs <- table(alleles) / length(alleles)
  
  # Return the minor allele frequency
  return(min(allele_freqs))
}

# Step D: Define a function to calculate missing data rate 
calculate_missing_rate <- function(genotypes) {
  return(mean(is.na(genotypes) | genotypes %in% c("N", "-", "NN", "NA", "")))
}

# Step E: Apply the filters 
maf_values <- apply(geno_data, 1, calculate_maf)
missing_rates <- apply(geno_data, 1, calculate_missing_rate)

# Combine results with original data
hapmap$MAF <- maf_values
hapmap$MissingRate <- missing_rates

# Filter for MAF >= 0.05 and missing rate <= 0.10
filtered_hapmap <- hapmap %>%
  filter(MAF >= 0.05, MissingRate <= 0.10)

# Step F: Save the filtered HapMap 
write.table(filtered_hapmap[, 1:(ncol(hapmap) - 2)], 
            file = "filtered_MAF_hapmap.hmp.txt", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE) 

