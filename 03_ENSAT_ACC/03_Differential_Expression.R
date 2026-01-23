# -------------------------------------------------------------------------
# Script: 03_Differential_Expression.R
# Description: Statistical comparison (C1A vs C1B) for HLA and Immune Cells.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, coin, openxlsx, here)

# Load Data
processed_path <- here("03_ENSAT_ACC", "processed_data")
load(file.path(processed_path, "ENSAT_HLA_Clinical.RData"))
gexp_hla <- as.data.frame(t(ENSAT[[1]])) # Transpose: Samples as rows
meta <- ENSAT[[2]]

mcp_scores <- read.xlsx(file.path(processed_path, "ENSAT_MCPcounter_scores.xlsx"))
rownames(mcp_scores) <- mcp_scores$ID

# 1. Prepare Stats Function
run_wilcox <- function(df, group_col, group1, group2) {
  df$Group <- meta[[group_col]][match(rownames(df), rownames(meta))]
  df <- df %>% filter(Group %in% c(group1, group2))
  
  res <- data.frame()
  features <- setdiff(colnames(df), "Group")
  
  for(f in features) {
    g1 <- df[df$Group == group1, f]
    g2 <- df[df$Group == group2, f]
    
    wt <- tryCatch(wilcox.test(g1, g2), error=function(e) list(p.value=NA))
    
    res <- rbind(res, data.frame(
      Feature = f,
      Mean_G1 = mean(g1, na.rm=T),
      Mean_G2 = mean(g2, na.rm=T),
      P_Value = wt$p.value
    ))
  }
  res$FDR <- p.adjust(res$P_Value, method="BH")
  return(res)
}

# 2. Run for HLA (C1A vs C1B)
# Note: ENSAT Metadata uses "C1A" and "C1B" in `Consensus.clustering.K2`
res_hla <- run_wilcox(gexp_hla, "Consensus.clustering.K2", "C1A", "C1B")
write.xlsx(res_hla, file.path(processed_path, "Stats_HLA_C1A_vs_C1B.xlsx"))

# 3. Run for Immune Cells (Selected 4 types)
cells_interest <- c("T_cells_MCPcounter", "CD8_T_cells_MCPcounter", 
                    "Cytotoxic_lymphocytes_MCPcounter", "B_lineage_MCPcounter")
mcp_subset <- mcp_scores[, cells_interest]

res_mcp <- run_wilcox(mcp_subset, "Consensus.clustering.K2", "C1A", "C1B")
write.xlsx(res_mcp, file.path(processed_path, "Stats_Immune_C1A_vs_C1B.xlsx"))

# 4. Generate Boxplots (Example for Immune)
pdf(file.path(processed_path, "Boxplots_Immune_C1A_C1B.pdf"))
mcp_subset$Group <- meta$Consensus.clustering.K2[match(rownames(mcp_subset), rownames(meta))]
mcp_long <- mcp_subset %>% pivot_longer(-Group, names_to="Cell", values_to="Score")

print(
  ggplot(mcp_long, aes(x=Group, y=Score, fill=Group)) +
    geom_boxplot(outlier.shape=NA) +
    geom_jitter(width=0.2, alpha=0.5) +
    facet_wrap(~Cell, scales="free") +
    scale_fill_manual(values=c("C1A"="#e41a1c", "C1B"="#377eb8")) +
    theme_bw() +
    labs(title="Immune Infiltration: C1A vs C1B")
)
dev.off()