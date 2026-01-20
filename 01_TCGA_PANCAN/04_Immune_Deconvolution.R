# -------------------------------------------------------------------------
# Project: Supplementary Material - Giner et al. (2025)
# Cohort: TCGA-PANCAN (ACC Subset)
# Script: 04_Immune_Deconvolution.R
# Description: Immune cell infiltration analysis using EPIC (IOBR package).
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, IOBR, ComplexHeatmap, circlize, RColorBrewer, here)

# Load Data
load(here("data", "coldataACC.RData"))
# Load the pre-extracted ACC Whole Transcriptome from Script 01
load(here("01_TCGA_PANCAN", "processed_data", "ACC_Whole_Transcriptome.RData"))

output_dir <- here("01_TCGA_PANCAN", "results", "immune_deconvolution")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Prepare Matrix
# rnaseq_acc_whole has 'gene_id' and samples.
# Need to convert Ensembl ID to Symbol for EPIC?
# Usually EPIC/IOBR works best with Symbols or TPM. 
# Assuming conversion logic is needed similar to Script 01 but for all genes.
# For simplicity, let's assume Ensembl IDs are handled or we load annotation.
# *Correction*: In original script, you loaded Annotation again. Let's do a quick mapping.

# Load Annotation (Lightweight part)
annotation <- rtracklayer::import(here("data", "gencode.v23.annotation.gff3.gz"))
gene_map <- as.data.frame(annotation) %>% 
  filter(type == "gene") %>%
  mutate(gene_id = str_remove(gene_id, "\\..*")) %>%
  select(gene_id, gene_name) %>%
  distinct()

# Map IDs to Symbols
eset_acc <- rnaseq_acc_whole %>%
  inner_join(gene_map, by = "gene_id") %>%
  select(-gene_id) %>%
  # Handle duplicates: keep max expression or mean
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), max)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

# -------------------------------------------------------------------------
# Run EPIC Deconvolution
# -------------------------------------------------------------------------
message("Running EPIC deconvolution...")
epic_res <- deconvo_tme(eset = eset_acc, method = "epic", arrays = FALSE)

# -------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------
# Align with Metadata
coldataACC$barcode <- substr(coldataACC$barcode, 1, 15)
epic_plot <- epic_res %>%
  filter(ID %in% coldataACC$barcode) %>%
  column_to_rownames("ID")

# Remove "other cells" if present and scale
if("otherCells" %in% colnames(epic_plot)) epic_plot <- epic_plot %>% select(-otherCells)
mat_epic <- t(scale(log2(epic_plot + 0.01)))

# Align order
common_samples <- intersect(colnames(mat_epic), coldataACC$barcode)
mat_epic <- mat_epic[, common_samples]
meta_ordered <- coldataACC[match(common_samples, coldataACC$barcode), ]

# Heatmap
pdf(file.path(output_dir, "EPIC_Infiltration_Heatmap.pdf"), width = 10, height = 6)
Heatmap(mat_epic,
        name = "Infiltration (Z-score)",
        column_split = meta_ordered$steroid,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        col = colorRamp2(seq(-2, 2, length=9), rev(brewer.pal(9, "RdYlBu"))))
dev.off()

message("Deconvolution complete.")