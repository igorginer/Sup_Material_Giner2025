# -------------------------------------------------------------------------
# Script: 07_Detailed_Heatmaps.R
# Description: Generates the comprehensive vertical heatmaps (HLA + T/B + Exhaustion)
#              stratified by Steroid Phenotype.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ComplexHeatmap, circlize, RColorBrewer, tidyverse, here)

# Load Processed Data
processed_path <- here("02_TCGA_ACC", "processed_data")
load(file.path(processed_path, "coldataACC_Clean.RData"))
load(file.path(processed_path, "rnaseq_hla.RData"))
load(file.path(processed_path, "rnaseq_exhaustion.RData"))
load(file.path(processed_path, "rnaseq_tb.RData"))
load(file.path(processed_path, "annotation_clean.RData"))

output_dir <- here("02_TCGA_ACC", "results", "complex_heatmaps")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Helper to process matrix
prep_mat <- function(df, annot_df, id_col="gene_id", name_col="gene_name", samples) {
  # Filter samples
  df_filt <- df %>% select(all_of(id_col), any_of(samples))
  
  # Map ID to Name if needed
  if(id_col == "gene_id") {
    df_filt <- df_filt %>% left_join(annot_df[,c(id_col, name_col)], by=id_col)
    mat <- df_filt %>% select(-all_of(id_col)) %>% column_to_rownames(name_col) %>% as.matrix()
  } else {
    mat <- df_filt %>% column_to_rownames(name_col) %>% as.matrix()
  }
  
  # Z-score & Clamp
  mat <- 2^mat - 1 # Unlog
  mat <- t(scale(t(mat)))
  mat[mat > 2] <- 2
  mat[mat < -2] <- -2
  return(mat)
}

# -------------------------------------------------------------------------
# Generate Plot for a Group (HSP or LSP)
# -------------------------------------------------------------------------
generate_complex_heatmap <- function(group_name) {
  
  # 1. Select Samples
  meta_sub <- coldataACC %>% filter(steroid == group_name)
  samps <- meta_sub$barcode
  # Ensure samples match column format (check if 12 or 15 chars in RNA)
  # Assuming RNA uses 15 chars based on prep script
  
  # 2. Prep Matrices
  mat_hla <- prep_mat(rnaseq_hla, annot_df, "gene_id", "gene_name", samps)
  mat_exh <- prep_mat(rnaseq_exhaustion, annot_df, "gene_id", "gene_name", samps)
  mat_tb  <- prep_mat(rnaseq_tb, annot_df, "gene_id", "gene_name", samps)
  
  # Align samples across matrices
  common <- intersect(colnames(mat_hla), colnames(mat_exh))
  common <- intersect(common, colnames(mat_tb))
  
  mat_hla <- mat_hla[, common]
  mat_exh <- mat_exh[, common]
  mat_tb  <- mat_tb[, common]
  meta_ord <- meta_sub[match(common, meta_sub$barcode), ]
  
  # 3. Annotation
  ha <- HeatmapAnnotation(
    Cortisol = meta_ord$cortisol.excess,
    Stage = meta_ord$tumor_stage,
    col = list(Cortisol = c("Cortisol"="#e41a1c", "No"="#377eb8"),
               Stage = c("stage i"="#bdbdbd", "stage ii"="#969696", "stage iii"="#353535", "stage iv"="#000000"))
  )
  
  # 4. Draw
  col_fun = colorRamp2(seq(-2, 2, length=11), rev(brewer.pal(11, "RdYlBu")))
  
  ht1 <- Heatmap(mat_hla, name="HLA", top_annotation = ha, col=col_fun, cluster_columns=TRUE, show_column_names=FALSE)
  ht2 <- Heatmap(mat_tb, name="T/B Markers", col=col_fun, cluster_columns=FALSE, show_column_names=FALSE)
  ht3 <- Heatmap(mat_exh, name="Exhaustion", col=col_fun, cluster_columns=FALSE, show_column_names=FALSE)
  
  pdf(file.path(output_dir, paste0("Complex_Heatmap_", group_name, ".pdf")), width=10, height=14)
  draw(ht1 %v% ht2 %v% ht3)
  dev.off()
}

generate_complex_heatmap("Steroid_High")
generate_complex_heatmap("Steroid_Low")