# -------------------------------------------------------------------------
# Project: Supplementary Material - Giner et al. (2025)
# Cohort: TCGA-ACC
# Script: 01_Data_Preparation.R
# Description: Imports raw counts, clinical metadata, and extracts specific
#              gene subsets (HLA, Exhaustion, T/B markers) for downstream analysis.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rtracklayer, stringr, here)

# Define paths
raw_data_path <- here("data")
processed_path <- here("02_TCGA_ACC", "processed_data")
if (!dir.exists(processed_path)) dir.create(processed_path, recursive = TRUE)

# -------------------------------------------------------------------------
# 1. Load Metadata
# -------------------------------------------------------------------------
message("Loading Metadata...")
load(file.path(raw_data_path, "coldataACC.RData"))
coldataACC$barcode <- substr(coldataACC$barcode, 1, 15)

# Filter samples with Steroid Info
coldataACC <- coldataACC[!is.na(coldataACC$steroid), ]

# -------------------------------------------------------------------------
# 2. Load Annotation
# -------------------------------------------------------------------------
# Note: Using gencode.v36 as per your original script
annotation <- import(file.path(raw_data_path, "gencode.v23.annotation.gff3.gz")) # Use the version available in your data folder
annot_df <- as.data.frame(annotation) %>% 
  filter(type == "gene") %>%
  mutate(gene_id = str_remove(gene_id, "\\..*"))

# -------------------------------------------------------------------------
# 3. Load Raw RNAseq (TCGA-ACC)
# -------------------------------------------------------------------------
message("Loading Expression Data...")
rnaseq <- read_tsv(file.path(raw_data_path, "TCGA-ACC.star_tpm.tsv.gz"))
colnames(rnaseq)[1] <- "gene_id"
rnaseq$gene_id <- str_remove(rnaseq$gene_id, "\\..*")

# -------------------------------------------------------------------------
# 4. Extract Gene Subsets
# -------------------------------------------------------------------------

# A. HLA and B2M Genes
hla_genes_df <- annot_df %>% 
  filter(str_detect(gene_name, "^HLA|^B2M")) %>%
  select(gene_name, gene_id) %>% distinct()

rnaseq_hla <- rnaseq %>% 
  filter(gene_id %in% hla_genes_df$gene_id) %>%
  left_join(hla_genes_df, by = "gene_id") %>%
  select(gene_name, gene_id, everything())

# B. Immune Markers (Exhaustion & T/B Cells)
# IDs sourced from your Script 21/22
exhaustion_ids <- c("ENSG00000113580","ENSG00000188389", "ENSG00000135077","ENSG00000089692",
                    "ENSG00000181847","ENSG00000136634", "ENSG00000138185","ENSG00000198846",
                    "ENSG00000165030","ENSG00000081059", "ENSG00000163599")

tb_ids <- c("ENSG00000198851", "ENSG00000160654", "ENSG00000167286", "ENSG00000188404", "ENSG00000126353", 
            "ENSG00000153563", "ENSG00000172116", "ENSG00000111537", "ENSG00000232810", "ENSG00000163508", 
            "ENSG00000010610", "ENSG00000186810", "ENSG00000109471", "ENSG00000113520", "ENSG00000168811", 
            "ENSG00000150782", "ENSG00000115415", "ENSG00000138378", "ENSG00000183813", "ENSG00000178562", 
            "ENSG00000081237", "ENSG00000073861", "ENSG00000121966", "ENSG00000113249", "ENSG00000107485", 
            "ENSG00000112115", "ENSG00000105329", "ENSG00000092969", "ENSG00000119699", "ENSG00000136244", 
            "ENSG00000177875", "ENSG00000138684", "ENSG00000134460", "ENSG00000168685", "ENSG00000049768", 
            "ENSG00000126561", "ENSG00000026508", "ENSG00000211899", "ENSG00000204287", "ENSG00000177455", 
            "ENSG00000156738", "ENSG00000101017", "ENSG00000120949", "ENSG00000004468", "ENSG00000115884", 
            "ENSG00000211896", "ENSG00000048462", "ENSG00000117322", "ENSG00000139193", "ENSG00000110848", 
            "ENSG00000026103", "ENSG00000158473", "ENSG00000197471", "ENSG00000104921", "ENSG00000211898", 
            "ENSG00000110448")

rnaseq_exhaustion <- rnaseq %>% filter(gene_id %in% exhaustion_ids)
rnaseq_tb <- rnaseq %>% filter(gene_id %in% tb_ids)

# -------------------------------------------------------------------------
# 5. Save Processed Data
# -------------------------------------------------------------------------
message("Saving RData objects...")

save(coldataACC, file = file.path(processed_path, "coldataACC_Clean.RData"))
save(rnaseq_hla, file = file.path(processed_path, "rnaseq_hla.RData"))
save(rnaseq_exhaustion, file = file.path(processed_path, "rnaseq_exhaustion.RData"))
save(rnaseq_tb, file = file.path(processed_path, "rnaseq_tb.RData"))
save(annotation, annot_df, file = file.path(processed_path, "annotation_clean.RData"))

message("Data preparation complete.")