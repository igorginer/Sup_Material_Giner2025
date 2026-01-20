# -------------------------------------------------------------------------
# Project: Supplementary Material - Giner et al. (2025)
# Cohort: TCGA-PANCAN
# Script: 02_Global_HLA_Analysis.R
# Description: Generates heatmaps of median expression and boxplots 
#              comparing HLA expression across cancer types.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, pheatmap, reshape2, here)

# Load Processed Data
load(here("01_TCGA_PANCAN", "processed_data", "PanCan_HLA_Medians.RData"))
load(here("01_TCGA_PANCAN", "processed_data", "PanCan_HLA_Expression.RData"))
load(here("01_TCGA_PANCAN", "processed_data", "Sample_Lists.RData"))

output_dir <- here("01_TCGA_PANCAN", "results", "global_analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Genes of Interest (Schaafsma et al., 2021)
hlas_interest <- c("HLA-DPA1", "HLA-DRB1", "HLA-DPB1", "HLA-DRA", "HLA-E", 
                   "HLA-B", "HLA-A", "HLA-C", "B2M", "HLA-DQA1", "HLA-DRB5", 
                   "HLA-DQB1", "HLA-DMA", "HLA-DMB", "HLA-F", "HLA-DRB6", 
                   "HLA-DQA2", "HLA-DOA", "HLA-DQB2", "HLA-DOB", "HLA-G", "HLA-DPB2")

# -------------------------------------------------------------------------
# Part 1: Heatmaps of Medians
# -------------------------------------------------------------------------
plot_median_heatmap <- function(data_df, filename, normalized = FALSE) {
  
  # Filter genes
  data_filt <- data_df %>% filter(gene_name %in% hlas_interest)
  mat <- as.matrix(data_filt[, 3:ncol(data_filt)])
  rownames(mat) <- data_filt$gene_name
  
  if (normalized) {
    # Z-score standardization
    mat <- t(scale(t(mat)))
    # Clamp values (-2 to 2)
    mat[mat > 2] <- 2
    mat[mat < -2] <- -2
  }
  
  pdf(file.path(output_dir, filename), width = 10, height = 8)
  pheatmap(mat, 
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           show_colnames = TRUE,
           main = ifelse(normalized, "Normalized Median Expression", "Raw Median Expression"))
  dev.off()
}

# Generate Plots
plot_median_heatmap(medians_general, "Heatmap_Median_General_Raw.pdf", normalized = FALSE)
plot_median_heatmap(medians_acc_split, "Heatmap_Median_ACC_Split_Raw.pdf", normalized = FALSE)
plot_median_heatmap(medians_general, "Heatmap_Median_General_Normalized.pdf", normalized = TRUE)
plot_median_heatmap(medians_acc_split, "Heatmap_Median_ACC_Split_Normalized.pdf", normalized = TRUE)

# -------------------------------------------------------------------------
# Part 2: Boxplots Comparison (Contextualizing ACC)
# -------------------------------------------------------------------------
message("Generating boxplots...")

# Filter Expression Data
rnaseq_filt <- rnaseq_hla %>% filter(gene_name %in% hlas_interest)

# Helper function to prepare long format data
prepare_plot_data <- function(sample_list, expr_data) {
  cancer_df <- map_dfr(names(sample_list), function(ctype) {
    tibble(sample = sample_list[[ctype]], cancer_type = ctype)
  })
  
  # Reshape expression data
  expr_long <- expr_data %>%
    pivot_longer(cols = -c(gene_name, gene_id), names_to = "sample", values_to = "expression") %>%
    inner_join(cancer_df, by = "sample")
  
  return(expr_long)
}

# Generate Data for ACC Split scenario
plot_data <- prepare_plot_data(tumor_samples_acc_split, rnaseq_filt)

# Loop through genes
for (gene in hlas_interest) {
  
  gene_df <- plot_data %>% filter(gene_name == gene)
  
  # Order cancer types by median expression
  medians <- gene_df %>% group_by(cancer_type) %>% summarize(med = median(expression, na.rm=T))
  gene_df$cancer_type <- factor(gene_df$cancer_type, levels = medians$cancer_type[order(medians$med)])
  
  # Plot
  p <- ggplot(gene_df, aes(x = cancer_type, y = expression, fill = cancer_type)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = c("ACC_High_Steroid" = "#F8766D", "ACC_Low_Steroid" = "#00BFC4"), 
                      na.value = "grey90") +
    labs(title = paste(gene, "- PanCancer Expression"), y = "TPM", x = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")
  
  ggsave(file.path(output_dir, paste0("Boxplot_", gene, "_PanCan.pdf")), p, width = 12, height = 6)
}

message("Visualization complete.")