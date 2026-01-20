# -------------------------------------------------------------------------
# Project: Supplementary Material - Giner et al. (2025)
# Cohort: TCGA-ACC
# Script: 02_HLA_Clustering.R
# Description: Unsupervised hierarchical clustering of HLA/B2M expression.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ComplexHeatmap, circlize, RColorBrewer, tidyverse, here)

# Load Data
processed_path <- here("02_TCGA_ACC", "processed_data")
load(file.path(processed_path, "coldataACC_Clean.RData"))
load(file.path(processed_path, "rnaseq_hla.RData"))
load(here("data", "tcgaACC_pre_processed.RData")) # For C1A/C1B

# 1. Prepare Metadata
# Match C1A/C1B from pre-processed data
tcgaACC$barcode <- substr(tcgaACC$barcode, 1, 15)
idx <- match(coldataACC$barcode, tcgaACC$barcode)
coldataACC$paper_C1A.C1B <- as.factor(as.character(tcgaACC$paper_C1A.C1B[idx]))

# Set final barcode length to 12 for RNAseq matching
coldataACC$barcode_12 <- substr(coldataACC$barcode, 1, 12)
rownames(coldataACC) <- coldataACC$barcode_12

# 2. Prepare Expression Matrix
hlas_of_interest <- c("HLA-DPA1", "HLA-DRB1", "HLA-DPB1", "HLA-DRA", "HLA-E", 
                      "HLA-B", "HLA-A", "HLA-C", "B2M", "HLA-DQA1", "HLA-DRB5", 
                      "HLA-DQB1", "HLA-DMA", "HLA-DMB", "HLA-F", "HLA-DRB6", 
                      "HLA-DQA2", "HLA-DOA", "HLA-DQB2", "HLA-DOB", "HLA-G", "HLA-DPB2")

mat <- rnaseq_hla %>%
  filter(gene_name %in% hlas_of_interest) %>%
  select(gene_name, any_of(coldataACC$barcode_12)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

# 3. Standardization (Z-score)
mat <- 2^mat - 1 # Unlog first (if source was log2)
mat_scaled <- t(scale(t(mat)))
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# 4. Annotation
common_samples <- intersect(colnames(mat_scaled), rownames(coldataACC))
coldata_ordered <- coldataACC[common_samples, ]
mat_scaled <- mat_scaled[, common_samples]

ha <- HeatmapAnnotation(
  steroid = coldata_ordered$steroid,
  C1A.C1B = coldata_ordered$paper_C1A.C1B,
  cortisol = coldata_ordered$cortisol.excess,
  stage = coldata_ordered$tumor_stage,
  col = list(
    steroid = c("Steroid_High"="#e41a1c", "Steroid_Low"="#377eb8"),
    C1A.C1B = c("C1A"="#e41a1c", "C1B"="#377eb8"),
    cortisol = c("Cortisol"="#e41a1c", "No"="#377eb8"),
    stage = c("stage i"= "#bdbdbd", "stage ii"="#969696", "stage iii"="#353535", "stage iv"="#000000")
  )
)

# 5. Plot
output_dir <- here("02_TCGA_ACC", "results", "clustering")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pdf(file.path(output_dir, "Global_HLA_Clustering.pdf"), width = 12, height = 9)
Heatmap(mat_scaled,
        name = "Z-score",
        top_annotation = ha,
        cluster_columns = TRUE,
        column_split = 4,
        cluster_rows = FALSE, # Keep custom order if desired, or TRUE
        show_column_names = FALSE,
        col = colorRamp2(seq(-2, 2, length=11), rev(brewer.pal(11, "RdYlBu"))))
dev.off()