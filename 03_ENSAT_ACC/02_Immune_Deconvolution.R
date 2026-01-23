# -------------------------------------------------------------------------
# Script: 02_Immune_Deconvolution.R
# Description: Prepares full transcriptome data for MCP-counter and runs it via IOBR.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, GEOquery, Biobase, hugene20sttranscriptcluster.db, AnnotationHub, IOBR, openxlsx, here)

# Paths
processed_path <- here("03_ENSAT_ACC", "processed_data")
if (!dir.exists(processed_path)) dir.create(processed_path)

# 1. Download Full Data Again (Need full matrix for deconvolution, not just HLA)
message("Fetching full expression matrix...")
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*4)
gds <- getGEO("GSE49278")
e <- gds[[1]]
gexp_full <- exprs(e)
pdat <- pData(e)

# 2. Get Protein Coding Genes
message("Fetching protein-coding gene annotations...")
ah <- AnnotationHub()
# Query specific version to ensure reproducibility (Homo sapiens EnsDb 103)
edb <- query(ah, pattern = c("Homo sapiens", "EnsDb", 103))[[1]]
gns <- genes(edb, return.type= "data.frame")
protein_coding_sym <- gns$symbol[gns$gene_biotype == "protein_coding"]

# 3. Map Probes
gene_names <- select(hugene20sttranscriptcluster.db, rownames(gexp_full), c("PROBEID","SYMBOL"))
gene_names <- gene_names %>% 
  filter(!is.na(SYMBOL)) %>%
  filter(SYMBOL %in% protein_coding_sym)

# 4. Collapse Probes (Max Mean method for simplicity in deconvolution prep)
# Using IOBR's built-in tool or simple aggregation
message("Collapsing probes to gene symbols...")
gexp_subset <- gexp_full[gene_names$PROBEID, ]
rownames(gexp_subset) <- gene_names$SYMBOL[match(rownames(gexp_subset), gene_names$PROBEID)]

# Handle duplicates (take max)
gexp_mcp <- gexp_subset %>% 
  as.data.frame() %>% 
  rownames_to_column("Gene") %>% 
  group_by(Gene) %>% 
  summarise(across(everything(), max)) %>% 
  column_to_rownames("Gene") %>% 
  as.matrix()

# 5. Run MCP-counter
message("Running MCP-counter...")
mcp_res <- deconvo_tme(eset = gexp_mcp, method = "mcpcounter")

# 6. Save Results
write.xlsx(mcp_res, file.path(processed_path, "ENSAT_MCPcounter_scores.xlsx"))
message("Deconvolution complete.")