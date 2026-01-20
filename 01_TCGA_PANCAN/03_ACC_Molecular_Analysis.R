# -------------------------------------------------------------------------
# Project: Supplementary Material - Giner et al. (2025)
# Cohort: TCGA-PANCAN (ACC Subset)
# Script: 03_ACC_Molecular_Analysis.R
# Description: Detailed clustering of ACC samples and statistical testing
#              (HSP vs LSP) for HLA genes.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ComplexHeatmap, circlize, RColorBrewer, coin, openxlsx, here)

# Load Data
load(here("data", "coldataACC.RData"))
load(here("data", "tcgaACC_pre_processed.RData"))
load(here("01_TCGA_PANCAN", "processed_data", "PanCan_HLA_Expression.RData"))

output_dir <- here("01_TCGA_PANCAN", "results", "acc_molecular")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Prepare ACC Data
coldataACC$barcode <- substr(coldataACC$barcode, 1, 15)
tcgaACC$barcode <- substr(tcgaACC$barcode, 1, 15)

# Subset Expression for ACC samples only
acc_samples <- coldataACC$barcode
rnaseq_acc <- rnaseq_hla %>% 
  select(gene_name, any_of(acc_samples)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

# Align Metadata
meta_acc <- coldataACC %>% 
  filter(barcode %in% colnames(rnaseq_acc)) %>%
  left_join(tcgaACC %>% select(barcode, paper_C1A.C1B, Leukocyte.Fraction), by = "barcode") %>%
  mutate(steroid = factor(steroid, levels = c("Steroid_Low", "Steroid_High"))) %>%
  column_to_rownames("barcode")

# Ensure order matches
rnaseq_acc <- rnaseq_acc[, rownames(meta_acc)]

# -------------------------------------------------------------------------
# Part 1: Complex Heatmap
# -------------------------------------------------------------------------
# Standardization (Z-score)
mat_scaled <- t(scale(t(rnaseq_acc)))
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# Annotations
ha <- HeatmapAnnotation(
  Steroid = meta_acc$steroid,
  Subtype = meta_acc$paper_C1A.C1B,
  Cortisol = meta_acc$cortisol.excess,
  Stage = meta_acc$tumor_stage,
  col = list(
    Steroid = c("Steroid_High"="#e41a1c", "Steroid_Low"="#377eb8"),
    Subtype = c("C1A"="#e41a1c", "C1B"="#377eb8"),
    Cortisol = c("Cortisol"="black", "No"="#525252")
  )
)

pdf(file.path(output_dir, "ACC_Detailed_Heatmap.pdf"), width = 10, height = 8)
Heatmap(mat_scaled, 
        name = "Z-score",
        top_annotation = ha,
        column_split = meta_acc$steroid,
        show_column_names = FALSE,
        col = colorRamp2(seq(-2, 2, length=11), rev(brewer.pal(11, "RdYlBu"))))
dev.off()

# -------------------------------------------------------------------------
# Part 2: Statistical Analysis (HSP vs LSP)
# -------------------------------------------------------------------------
stats_results <- data.frame()

# Define groups
hsp_ids <- rownames(meta_acc)[meta_acc$steroid == "Steroid_High"]
lsp_ids <- rownames(meta_acc)[meta_acc$steroid == "Steroid_Low"]

genes <- rownames(rnaseq_acc)

for (g in genes) {
  val_hsp <- rnaseq_acc[g, hsp_ids]
  val_lsp <- rnaseq_acc[g, lsp_ids]
  
  # Skip if all zero
  if(sum(val_hsp) == 0 && sum(val_lsp) == 0) next
  
  # Tests
  t_res <- tryCatch(t.test(val_hsp, val_lsp), error=function(e) list(p.value=NA))
  w_res <- tryCatch(wilcox.test(val_hsp, val_lsp), error=function(e) list(p.value=NA))
  
  # Permutation (Coin)
  df_perm <- data.frame(expr = c(val_hsp, val_lsp), 
                        grp = factor(c(rep("HSP", length(val_hsp)), rep("LSP", length(val_lsp)))))
  perm_res <- tryCatch(pvalue(oneway_test(expr ~ grp, data=df_perm, distribution="approximate")), 
                       error=function(e) NA)
  
  # Save
  stats_results <- rbind(stats_results, data.frame(
    gene = g,
    mean_HSP = mean(val_hsp, na.rm=T),
    mean_LSP = mean(val_lsp, na.rm=T),
    p_ttest = t_res$p.value,
    p_wilcox = w_res$p.value,
    p_perm = perm_res
  ))
  
  # Plot Boxplot (Optional, saves space not to plot all 22 here, but can enable)
}

write.csv(stats_results, file.path(output_dir, "ACC_HLA_Stats.csv"), row.names = FALSE)
message("Analysis complete.")