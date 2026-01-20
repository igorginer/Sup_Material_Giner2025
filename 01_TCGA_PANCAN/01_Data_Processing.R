# -------------------------------------------------------------------------
# Project: Supplementary Material - Giner et al. (2025)
# Cohort: TCGA-PANCAN
# Script: 01_Data_Processing.R
# Description: Imports raw clinical and RNAseq data, filters tumor samples,
#              extracts HLA/B2M genes, and saves processed RData objects.
# -------------------------------------------------------------------------

# 1. Setup
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, vroom, rtracklayer, stringr, here)

# Define paths relative to the project root
# Ensure raw data is in: project_root/data/
raw_data_path <- here("data")
processed_data_path <- here("01_TCGA_PANCAN", "processed_data")
if (!dir.exists(processed_data_path)) dir.create(processed_data_path)

# -------------------------------------------------------------------------
# 2. Clinical Data Processing
# -------------------------------------------------------------------------
message("Loading and processing clinical data...")

# Load PanCancer Survival Data (Clinical)
clinical_general <- read_tsv(file.path(raw_data_path, "Survival_SupplementalTable_S1_20171025_xena_sp"))

# Load ACC Metadata (Local cohort)
load(file.path(raw_data_path, "coldataACC.RData"))
coldataACC$barcode <- substr(coldataACC$barcode, 1, 15)

# Split ACC samples by Steroid Phenotype
acc_high <- coldataACC$barcode[coldataACC$steroid == "Steroid_High" & !is.na(coldataACC$steroid)]
acc_low  <- coldataACC$barcode[coldataACC$steroid == "Steroid_Low"  & !is.na(coldataACC$steroid)]

# Create Sample Lists per Cancer Type
# List 1: General (ACC as one group)
samples_general <- split(clinical_general$sample, clinical_general$`cancer type abbreviation`)

# List 2: Specific (ACC split by steroid)
samples_acc_split <- samples_general
samples_acc_split[["ACC"]] <- NULL # Remove generic ACC
samples_acc_split[["ACC_High_Steroid"]] <- acc_high
samples_acc_split[["ACC_Low_Steroid"]]  <- acc_low
samples_acc_split <- samples_acc_split[sort(names(samples_acc_split))]

# Function to filter Tumor Samples (01-09)
get_tumor_samples <- function(sample_list) {
  lapply(sample_list, function(x) grep("-0[1-9]$", x, value = TRUE))
}

tumor_samples_general <- get_tumor_samples(samples_general)
tumor_samples_acc_split <- get_tumor_samples(samples_acc_split)

# -------------------------------------------------------------------------
# 3. RNAseq Data Processing (Heavy Step)
# -------------------------------------------------------------------------
message("Loading massive RNAseq file (this may take a while)...")

# Load Annotation first to identify HLA genes
annotation <- import(file.path(raw_data_path, "gencode.v23.annotation.gff3.gz"))
annot_df <- as.data.frame(annotation) %>% 
  filter(type == "gene") %>%
  mutate(gene_id = str_remove(gene_id, "\\..*"))

# Define HLA/B2M genes of interest
hla_genes <- annot_df %>% 
  filter(str_detect(gene_name, "^HLA|^B2M")) %>%
  select(gene_name, gene_id) %>% 
  distinct()

# Read RNAseq (Chunked reading not necessary if RAM > 16GB, but used vroom for speed)
# IMPORTANT: Ensure the file name matches your local file
rnaseq_raw <- vroom(file.path(raw_data_path, "tcga_RSEM_gene_tpm.gz"), 
                    col_types = cols(), show_col_types = FALSE)
colnames(rnaseq_raw)[1] <- "gene_id"
rnaseq_raw$gene_id <- str_remove(rnaseq_raw$gene_id, "\\..*")

# -------------------------------------------------------------------------
# 4. Extract Subsets and Save
# -------------------------------------------------------------------------

# A. Extract HLA Genes (Pan-Cancer)
message("Extracting HLA subset...")
rnaseq_hla <- rnaseq_raw %>% 
  filter(gene_id %in% hla_genes$gene_id) %>%
  left_join(hla_genes, by = "gene_id") %>%
  select(gene_name, gene_id, everything())

# B. Extract ACC Whole Transcriptome (For Deconvolution/Script 04)
# This avoids reloading the big file later
message("Extracting ACC Whole Transcriptome subset...")
acc_all_samples <- clinical_general$sample[clinical_general$`cancer type abbreviation` == "ACC"]
acc_samples_in_rna <- intersect(acc_all_samples, colnames(rnaseq_raw))
rnaseq_acc_whole <- rnaseq_raw %>% 
  select(gene_id, all_of(acc_samples_in_rna))

# C. Calculate Medians (Tumor Only)
message("Calculating medians per cancer type...")

calc_medians <- function(expr_data, sample_list) {
  # Initialize DF
  medians_df <- expr_data[, 1:2] # gene_name, gene_id
  
  for (cancer in names(sample_list)) {
    samps <- intersect(sample_list[[cancer]], colnames(expr_data))
    if (length(samps) > 0) {
      medians_df[[cancer]] <- apply(expr_data[, samps], 1, median, na.rm = TRUE)
    }
  }
  return(medians_df)
}

medians_general <- calc_medians(rnaseq_hla, tumor_samples_general)
medians_acc_split <- calc_medians(rnaseq_hla, tumor_samples_acc_split)

# -------------------------------------------------------------------------
# 5. Save Outputs
# -------------------------------------------------------------------------
message("Saving processed files...")

# Save Lists
save(tumor_samples_general, tumor_samples_acc_split, 
     file = file.path(processed_data_path, "Sample_Lists.RData"))

# Save Expression Data
save(rnaseq_hla, file = file.path(processed_data_path, "PanCan_HLA_Expression.RData"))
save(rnaseq_acc_whole, file = file.path(processed_data_path, "ACC_Whole_Transcriptome.RData"))

# Save Medians
save(medians_general, medians_acc_split, 
     file = file.path(processed_data_path, "PanCan_HLA_Medians.RData"))

message("Processing complete! Proceed to Script 02.")