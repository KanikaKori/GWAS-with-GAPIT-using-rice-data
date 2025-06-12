# This always gives you the latest version without needing devtools or Rtools
source("https://zzlab.net/GAPIT/gapit_functions.txt")
source("https://zzlab.net/GAPIT/emma.txt")

##If the Genotype IDs present in the Phenotype file and Genotype file is not same
# 1. Load your raw data
# Adjust file paths as needed
myY_raw <- read.table("Path_to_phenotype_file/pheno_data_for_gwas.txt", header = TRUE, sep = "\t") # Assuming tab-delimited
myGD_raw <- read_tsv("Path_to_genotype_file/Geno_data.hmp.txt", comment ="") # Assuming tab-delimited

# 2. Extract Taxa IDs from both files
# For phenotype file, it's the first column
pheno_ids <- as.character(myY_raw[, 1]) # Convert to character to be safe

# For genotype file (HapMap format), IDs are column names AFTER the first 11
# Make sure to handle the column names correctly based on your file
genotype_ids <- colnames(myGD_raw)[12:ncol(myGD_raw)] # Assuming 11 metadata columns

# 3. Find common IDs between the two sets
common_ids <- intersect(pheno_ids, genotype_ids)

# 4. Check how many common IDs you have
num_common_ids <- length(common_ids)
print(paste("Number of common individuals found:", num_common_ids))

# More robust way using the first column regardless of its name
myY_filtered <- myY_raw[myY_raw[[1]] %in% common_ids, ]

# 5. Subset the genotype data
# Select the first 11 metadata columns AND the columns corresponding to common_ids
# Need to select columns by name
myGD_filtered <- myGD_raw[, c(colnames(myGD_raw)[1:11], common_ids)]

# Reorder phenotype data for easy debugging
myY_filtered <- myY_filtered[match(common_ids, myY_filtered[[1]]), ]

# Reorder genotype data columns (only the taxa columns, keep metadata fixed)
myGD_filtered_taxa_only <- myGD_filtered[, common_ids]
myGD_filtered <- cbind(myGD_filtered[, 1:11], myGD_filtered_taxa_only)


# 6. Saving filtered phenotype data into a text file 
write.table(
  myY_filtered,
  file = "pheno_filtered.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 7. Saving filtered genotype data into a text file 
write.table(
  myGD_filtered,
  file = "geno_filtered.hmp.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 8.Run GAPIT with new filtered genotype and phenotype files
#Import data for gapit
myY <- read.table("pheno_filtered.txt", head = TRUE)
myG <- read.delim("geno_filtered.hmp.txt", head = FALSE)

#GWAS
myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3,
  model=c("GLM", "MLM", "FARMCPU", "CMLM"))
