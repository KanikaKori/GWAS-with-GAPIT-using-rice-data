# GWAS-with-GAPIT-using-rice-data
This R script performs a Genome-Wide Association Study (GWAS) using the GAPIT (Genome Association and Prediction Integrated Tool) pipeline. The workflow includes preprocessing, data filtering, and running multiple GWAS models using rice genotype and phenotype data.

## Steps and Explanation:

### 1. Load GAPIT functions
- GAPIT and EMMA scripts are sourced directly from the developers' server:

   source("https://zzlab.net/GAPIT/gapit_functions.txt")

   source("https://zzlab.net/GAPIT/emma.txt")

### 2. Load raw phenotype and genotype data
- The phenotype file is a tab-separated file containing traits with sample IDs.
- The genotype file is in HapMap format (.hmp.txt) with metadata in the first 11 columns and sample genotype data from column 12 onward.

### 3. Match sample IDs between phenotype and genotype
- The script extracts sample (taxa) IDs from both files and finds their intersection.
- Only the common individuals are kept for downstream analysis.

### 4. Filter the datasets
- Both phenotype and genotype files are filtered to include only the common IDs.
- Columns are reordered to match between files to avoid misalignment during GWAS.

### 5. Save the filtered data
- Two new files (pheno_filtered.txt and geno_filtered.hmp.txt) are saved and used as input for GAPIT.

### 6. Run GWAS using GAPIT
- The filtered data is loaded and GAPIT is run using multiple models:
   GLM – General Linear Model
   MLM – Mixed Linear Model
   FARMCPU – Fixed and random model circulating probability unification
   CMLM – Compressed MLM

### --Example command:

myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  PCA.total = 3,
  model = c("GLM", "MLM", "FARMCPU", "CMLM")
)

### 📁 Input File Format:
- Phenotype File: tab-separated, with the first column as genotype IDs (Taxa) and remaining columns as trait values.
- Genotype File: HapMap format with 11 metadata columns followed by genotype data for each Taxa.

### 📤 Output:
- GAPIT produces result files in the working directory, including:
  
   Manhattan plots
  
   QQ plots
  
   Summary tables of significant marker-trait associations for each model

