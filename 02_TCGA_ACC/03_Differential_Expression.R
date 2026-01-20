# -------------------------------------------------------------------------
# Script: 03_Differential_Expression.R
# Description: Wilcoxon/T-tests for HLA and Immune Markers across groups.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, coin, here, openxlsx)

load(here("02_TCGA_ACC", "processed_data", "rnaseq_hla.RData"))
load(here("02_TCGA_ACC", "processed_data", "rnaseq_exhaustion.RData"))
load(here("02_TCGA_ACC", "processed_data", "coldataACC_Clean.RData"))
load(here("data", "tcgaACC_pre_processed.RData")) # For C1A/C1B

output_dir <- here("02_TCGA_ACC", "results", "stats")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Helper Function for Stats
run_stats <- function(expr_df, metadata_df, group_col, group1, group2, id_col="gene_name") {
  
  # Format Expression
  expr_mat <- expr_df %>%
    select(all_of(id_col), any_of(metadata_df$barcode)) %>%
    column_to_rownames(id_col) %>%
    as.matrix() %>% t() %>% as.data.frame()
  
  # Add Group
  expr_mat$Group <- metadata_df[[group_col]][match(rownames(expr_mat), metadata_df$barcode)]
  expr_mat <- expr_mat %>% filter(Group %in% c(group1, group2))
  
  results <- data.frame()
  genes <- setdiff(colnames(expr_mat), "Group")
  
  for(g in genes) {
    g1_vals <- expr_mat[expr_mat$Group == group1, g]
    g2_vals <- expr_mat[expr_mat$Group == group2, g]
    
    # Wilcoxon
    w_test <- tryCatch(wilcox.test(g1_vals, g2_vals), error=function(e) list(p.value=NA))
    # T-test
    t_test <- tryCatch(t.test(g1_vals, g2_vals), error=function(e) list(p.value=NA))
    
    results <- rbind(results, data.frame(
      Gene = g,
      Mean_G1 = mean(g1_vals, na.rm=T),
      Mean_G2 = mean(g2_vals, na.rm=T),
      P_Wilcox = w_test$p.value,
      P_Ttest = t_test$p.value
    ))
  }
  
  results$FDR_Wilcox <- p.adjust(results$P_Wilcox, method = "BH")
  return(results)
}

# 1. HSP vs LSP (HLA)
coldataACC$barcode <- substr(coldataACC$barcode, 1, 15)
# Ensure RNAseq columns match metadata barcode length (likely 12 vs 15 check)
colnames(rnaseq_hla)[-c(1,2)] <- substr(colnames(rnaseq_hla)[-c(1,2)], 1, 15)

res_hla_steroid <- run_stats(rnaseq_hla, coldataACC, "steroid", "Steroid_High", "Steroid_Low")
write.csv(res_hla_steroid, file.path(output_dir, "HLA_Stats_HSP_vs_LSP.csv"))

# 2. C1A vs C1B (HLA)
tcgaACC$barcode <- substr(tcgaACC$barcode, 1, 15)
meta_merged <- coldataACC %>%
  left_join(as.data.frame(colData(tcgaACC)) %>% select(patient, paper_C1A.C1B), 
            by = c("barcode" = "patient"))

res_hla_subtype <- run_stats(rnaseq_hla, meta_merged, "paper_C1A.C1B", "C1A", "C1B")
write.csv(res_hla_subtype, file.path(output_dir, "HLA_Stats_C1A_vs_C1B.csv"))

# 3. Immune Markers (HSP vs LSP)
# Need to Map IDs to Symbols for Exhaustion first
load(here("02_TCGA_ACC", "processed_data", "annotation_clean.RData"))
rnaseq_exh_named <- rnaseq_exhaustion %>%
  left_join(annot_df %>% select(gene_id, gene_name), by="gene_id") %>%
  select(gene_name, everything(), -gene_id)

colnames(rnaseq_exh_named)[-1] <- substr(colnames(rnaseq_exh_named)[-1], 1, 15)
res_exh_steroid <- run_stats(rnaseq_exh_named, coldataACC, "steroid", "Steroid_High", "Steroid_Low")
write.csv(res_exh_steroid, file.path(output_dir, "Exhaustion_Stats_HSP_vs_LSP.csv"))