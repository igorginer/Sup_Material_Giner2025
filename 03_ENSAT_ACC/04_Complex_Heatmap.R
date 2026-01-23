# -------------------------------------------------------------------------
# Script: 04_Complex_Heatmap.R
# Description: Final visualization of C1A/C1B with Expression and MCP scores.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ComplexHeatmap, circlize, RColorBrewer, openxlsx, tidyverse, here)

# Load Data
processed_path <- here("03_ENSAT_ACC", "processed_data")
load(file.path(processed_path, "ENSAT_HLA_Clinical.RData"))
gexp_hla <- ENSAT[[1]] # Genes x Samples
meta <- ENSAT[[2]]
mcp <- read.xlsx(file.path(processed_path, "ENSAT_MCPcounter_scores.xlsx"), rowNames=TRUE)

# 1. Prepare Expression (Z-score)
mat_scaled <- t(scale(t(gexp_hla)))
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# 2. Prepare Annotation (Merge Meta + MCP)
mcp_subset <- mcp %>% select(T_cells_MCPcounter, CD8_T_cells_MCPcounter, 
                             Cytotoxic_lymphocytes_MCPcounter, B_lineage_MCPcounter)
# Scale MCP
mcp_scaled <- as.data.frame(scale(mcp_subset))
mcp_scaled[mcp_scaled > 2] <- 2
mcp_scaled[mcp_scaled < -2] <- -2

annot_df <- cbind(meta[, c("Consensus.clustering.K2", "OS")], mcp_scaled[rownames(meta), ])

# 3. Create Heatmap Function
make_heatmap <- function(group) {
  samps <- rownames(annot_df)[annot_df$Consensus.clustering.K2 == group]
  mat_sub <- mat_scaled[, samps]
  ann_sub <- annot_df[samps, ]
  
  ha <- HeatmapAnnotation(
    OS = as.factor(ann_sub$OS),
    T_Cells = ann_sub$T_cells_MCPcounter,
    CD8 = ann_sub$CD8_T_cells_MCPcounter,
    Cytotoxic = ann_sub$Cytotoxic_lymphocytes_MCPcounter,
    B_Cells = ann_sub$B_lineage_MCPcounter,
    col = list(
      OS = c("0"="#e6e6e6", "1"="black"),
      T_Cells = colorRamp2(c(-2, 2), c("white", "darkgreen")),
      CD8 = colorRamp2(c(-2, 2), c("white", "darkgreen")),
      Cytotoxic = colorRamp2(c(-2, 2), c("white", "darkgreen")),
      B_Cells = colorRamp2(c(-2, 2), c("white", "darkgreen"))
    )
  )
  
  Heatmap(mat_sub, name="Exp", top_annotation=ha, show_column_names=FALSE,
          col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
          column_title = paste("ENSAT -", group))
}

# 4. Save
pdf(here("03_ENSAT_ACC", "results", "Complex_Heatmaps.pdf"), width=10, height=8)
draw(make_heatmap("C1A"))
draw(make_heatmap("C1B"))
dev.off()