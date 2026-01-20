# -------------------------------------------------------------------------
# Script: 05_Immune_Correlations.R
# Description: Spearman correlations between HLA and Immune Markers/Scores.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(corrplot, RColorBrewer, tidyverse, here)

# Load Data
load(here("02_TCGA_ACC", "processed_data", "rnaseq_hla.RData"))
load(here("02_TCGA_ACC", "processed_data", "rnaseq_exhaustion.RData"))
load(here("02_TCGA_ACC", "processed_data", "coldataACC_Clean.RData"))

# Load External Data (Ensure these are in data folder)
mcp <- read.csv(here("data", "mcpcounter.csv"), row.names = 1) 
# Note: User provided 'mcpcounter.csv' in script 8. Assuming it exists.

output_dir <- here("02_TCGA_ACC", "results", "correlations")
if (!dir.exists(output_dir)) dir.create(output_dir)

# -------------------------------------------------------------------------
# 1. Prepare Matrices
# -------------------------------------------------------------------------
# Align samples
coldataACC$barcode <- substr(coldataACC$barcode, 1, 12) # Short barcode for correlation usually
rownames(mcp) <- substr(rownames(mcp), 1, 12)
rownames(mcp) <- gsub("\\.", "-", rownames(mcp))

# HLA Matrix
colnames(rnaseq_hla)[-c(1,2)] <- substr(colnames(rnaseq_hla)[-c(1,2)], 1, 12)
mat_hla <- rnaseq_hla %>% 
  select(-gene_id) %>% 
  column_to_rownames("gene_name") %>% 
  select(any_of(rownames(mcp))) %>% 
  as.matrix() %>% t()

# Filter Metadata
meta <- coldataACC %>% filter(barcode %in% rownames(mat_hla))

# -------------------------------------------------------------------------
# 2. Correlation Function (Robust)
# -------------------------------------------------------------------------
calc_corr_group <- function(hla_mat, immune_mat, group_samples, method="spearman") {
  
  # Subset samples
  samps <- intersect(group_samples, rownames(hla_mat))
  samps <- intersect(samps, rownames(immune_mat))
  
  m1 <- hla_mat[samps, ]
  m2 <- immune_mat[samps, ]
  
  # Cor test
  p_mat <- matrix(NA, ncol(m1), ncol(m2))
  r_mat <- matrix(NA, ncol(m1), ncol(m2))
  rownames(p_mat) <- rownames(r_mat) <- colnames(m1)
  colnames(p_mat) <- colnames(r_mat) <- colnames(m2)
  
  for(i in 1:ncol(m1)) {
    for(j in 1:ncol(m2)) {
      res <- cor.test(m1[,i], m2[,j], method=method)
      r_mat[i,j] <- res$estimate
      p_mat[i,j] <- res$p.value
    }
  }
  
  # Adjust P (Global)
  p_adj <- matrix(p.adjust(as.vector(p_mat), method="BH"), nrow=nrow(p_mat))
  return(list(r=r_mat, p_adj=p_adj))
}

# -------------------------------------------------------------------------
# 3. Run for HSP and LSP
# -------------------------------------------------------------------------
hsp_ids <- meta$barcode[meta$steroid == "Steroid_High"]
lsp_ids <- meta$barcode[meta$steroid == "Steroid_Low"]

# HLA vs MCP-Counter
res_hsp <- calc_corr_group(mat_hla, mcp, hsp_ids)
res_lsp <- calc_corr_group(mat_hla, mcp, lsp_ids)

# Plot
pdf(file.path(output_dir, "Correlation_HLA_MCP.pdf"), width = 10, height = 10)
corrplot(res_hsp$r, p.mat = res_hsp$p_adj, title = "HSP: HLA vs MCP", mar=c(0,0,2,0))
corrplot(res_lsp$r, p.mat = res_lsp$p_adj, title = "LSP: HLA vs MCP", mar=c(0,0,2,0))
dev.off()

# Save tables
write.csv(res_hsp$r, file.path(output_dir, "Corr_Rho_HSP.csv"))
write.csv(res_lsp$r, file.path(output_dir, "Corr_Rho_LSP.csv"))