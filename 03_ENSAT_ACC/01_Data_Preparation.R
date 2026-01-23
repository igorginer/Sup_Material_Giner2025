# -------------------------------------------------------------------------
# Project: Supplementary Material - Giner et al. (2025)
# Cohort: ENSAT (Validation)
# Script: 01_Data_Preparation.R
# Description: Downloads microarray data from GEO (GSE49278), maps probes 
#              to genes, constructs clinical metadata, and saves the object.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, GEOquery, Biobase, hugene20sttranscriptcluster.db, here)

# Define output path
processed_path <- here("03_ENSAT_ACC", "processed_data")
if (!dir.exists(processed_path)) dir.create(processed_path, recursive = TRUE)

# -------------------------------------------------------------------------
# 1. Download GEO Data (GSE49278)
# -------------------------------------------------------------------------
message("Downloading data from GEO (this may take a few minutes)...")
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*4)

# Download and extract expression set
gds <- getGEO("GSE49278")
e <- gds[[1]]
ENSAT_pdat <- pData(e)
gexp_ENSAT <- exprs(e)

# -------------------------------------------------------------------------
# 2. Probe to Gene Mapping & Cleaning
# -------------------------------------------------------------------------
message("Mapping probes to genes...")

gene.names <- select(hugene20sttranscriptcluster.db,
                     rownames(gexp_ENSAT), 
                     c("ENTREZID","SYMBOL","GENENAME"))
gene.names <- gene.names[complete.cases(gene.names$SYMBOL),]

# Filter for HLA and B2M only
rowdat.up_ENSAT <- gene.names %>%
  filter(grepl("^HLA-", SYMBOL) | SYMBOL == "B2M")

# Remove non-specific probes (mapped to multiple IDs)
duplicate_indices <- duplicated(rowdat.up_ENSAT$PROBEID) | duplicated(rowdat.up_ENSAT$PROBEID, fromLast = TRUE)
probes_to_remove <- unique(rowdat.up_ENSAT$PROBEID[duplicate_indices])

if (length(probes_to_remove) > 0) {
  rowdat.up_ENSAT <- rowdat.up_ENSAT[!rowdat.up_ENSAT$PROBEID %in% probes_to_remove, ]
}

# Resolve duplicated genes (multiple probes for one gene) by Max CV
gexp_subset <- as.data.frame(gexp_ENSAT[rowdat.up_ENSAT$PROBEID, ])
stats <- data.frame(
  PROBEID = rownames(gexp_subset),
  mean = apply(gexp_subset, 1, mean),
  sd = apply(gexp_subset, 1, sd)
)
stats$cv <- stats$sd / stats$mean

rowdat.up_ENSAT <- merge(rowdat.up_ENSAT, stats, by="PROBEID")

# Select probe with highest CV for each gene
rowdat.up_ENSAT <- rowdat.up_ENSAT %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = cv, n = 1) %>%
  ungroup()

# Create Final Expression Matrix
gexp_final <- gexp_ENSAT[rowdat.up_ENSAT$PROBEID, ]
rownames(gexp_final) <- rowdat.up_ENSAT$SYMBOL

# Filter for specific genes of interest (from your list)
genes_of_interest <- c("HLA-DPA1", "HLA-DRB1", "HLA-DPB1", "HLA-DRA", "HLA-E", 
                       "HLA-B", "HLA-A", "HLA-C", "B2M", "HLA-DQA1", "HLA-DRB5", 
                       "HLA-DQB1", "HLA-DMA", "HLA-DMB", "HLA-F", "HLA-DRB6", 
                       "HLA-DQA2", "HLA-DOA", "HLA-DQB2", "HLA-DOB", "HLA-G", "HLA-DPB2")

gexp_final <- gexp_final[rownames(gexp_final) %in% genes_of_interest, ]

# -------------------------------------------------------------------------
# 3. Construct Clinical Metadata (Hardcoded from Assié et al., 2014)
# -------------------------------------------------------------------------
message("Constructing clinical metadata...")

ENSAT_clinicaldat <- data.frame(
  stringsAsFactors = FALSE, check.names = FALSE,
  Sample = c("ACC1","ACC2","ACC3","ACC4","ACC5","ACC6","ACC7","ACC8","ACC9","ACC10","ACC11",
             "ACC12","ACC13","ACC14","ACC15","ACC16","ACC17","ACC18","ACC19","ACC20","ACC21",
             "ACC22","ACC23","ACC24","ACC25","ACC26","ACC27","ACC28","ACC29","ACC30","ACC31",
             "ACC32","ACC33","ACC34","ACC35","ACC36","ACC37","ACC38","ACC39","ACC40","ACC41",
             "ACC42","ACC43","ACC44","ACC45","ACC46","ACC47","ACC48","ACC49","ACC50","ACC51",
             "ACC52","ACC55"),
  Sex = c("F","F","M","F","M","F","F","F","M","F","M","F","F","F","F","F","F","F","F",
          "M","F","F","F","F","F","F","F","F","F","F","M","F","F","F","F","M","F","F","F","F",
          "F","M","M","F","F","F","F","M","F","F","M","F","M"),
  Age = c(70.3,25.6,40,53.3,72.9,18.3,77.5,50.7,63.9,27,29.6,79.3,46.2,43,53.9,45,41,37.2,
          81.6,67.5,42.3,39.7,25.2,41.7,37.9,23.9,59.5,75.5,37.6,34.1,26.1,26.5,48.4,58.8,
          49.6,54.3,79.6,29,44.5,28.5,68.9,28.9,52.4,30,46.3,59.4,18.6,39.6,40.2,53.8,30.6,
          44.6,52.9),
  Specific.Death = c("no","no","yes","no","yes","no","yes","yes","yes","no","yes","no","yes","yes",
                     "no","yes","no","yes","no","yes","no","no","yes","no","no","yes","no","no",
                     "no","yes","yes","yes","yes","yes","no","yes","yes","no","no","no","yes",
                     "yes","yes","no","no","no","yes","yes","no","no","no","no","yes"),
  Follow.up.months = c(151.8,131.1,23,147.8,0.5,142.8,5.2,36,0.6,73.6,74.4,12,12.7,9.5,115.1,85.1,
                       59.6,11.3,40.2,40.2,41.6,57.2,48.3,35.9,99.4,9.4,40.9,28,81.8,11.7,30.2,2,
                       19.8,34.3,129.8,57.9,7.6,118,119.5,144.3,19,21,24.3,111.5,97.6,12.9,55.4,
                       26,158.8,154.2,85.1,351.7,11.4),
  ENSAT.staging = c(2,1,4,2,4,1,2,2,4,2,4,2,4,3,2,NA,1,4,2,1,2,2,4,2,2,4,2,2,2,3,4,4,4,2,2,3,
                    4,2,2,2,3,4,4,2,2,2,2,2,2,4,2,NA,3)
)

ENSAT_c1a.c1b <- data.frame(
  stringsAsFactors = FALSE, check.names = FALSE,
  Sample = c("ACC1","ACC10","ACC11","ACC12", "ACC13","ACC14","ACC15","ACC16","ACC17","ACC18","ACC19",
             "ACC2","ACC20","ACC21","ACC22","ACC23","ACC24","ACC25","ACC26","ACC27","ACC28","ACC29",
             "ACC31","ACC32","ACC33","ACC35","ACC36","ACC37","ACC38","ACC39","ACC4","ACC40","ACC42",
             "ACC43","ACC44","ACC45","ACC46","ACC47","ACC48","ACC49","ACC5","ACC50","ACC52","ACC55",
             "ACC6","ACC8","ACC9"),
  Consensus.clustering.K2 = c("C1B","C1B","C1B","C1A","C1A", "C1A","C1B","C1A","C1B","C1A","C1A",
                              "C1B","C1A","C1B","C1B","C1A","C1A","C1B","C1A","C1A","C1B","C1B",
                              "C1A","C1A","C1A","C1B","C1A","C1A","C1A","C1A","C1B","C1B","C1A",
                              "C1A","C1A","C1B","C1A","C1A","C1A","C1B","C1A","C1B","C1B","C1A",
                              "C1A","C1A","C1B")
)

# Merge Clinical Data
rownames(ENSAT_clinicaldat) <- ENSAT_clinicaldat$Sample
rownames(ENSAT_c1a.c1b) <- ENSAT_c1a.c1b$Sample

# Filter to match GEO samples
samples_in_geo <- ENSAT_pdat$title
ENSAT_clinicaldat <- ENSAT_clinicaldat[samples_in_geo, ]
ENSAT_c1a.c1b <- ENSAT_c1a.c1b[samples_in_geo, ]

# Combine into final metadata
ENSAT_coldat <- cbind(ENSAT_pdat[,c(1,2)], ENSAT_clinicaldat, ENSAT_c1a.c1b)
ENSAT_coldat$OS <- ifelse(ENSAT_coldat$Specific.Death == "yes", 1, 0)
ENSAT_coldat$OS.time <- as.numeric(ENSAT_coldat$Follow.up.months)

# -------------------------------------------------------------------------
# 4. Save Outputs
# -------------------------------------------------------------------------
# Creating the List object as per your original logic
ENSAT <- list(gexp_final, ENSAT_coldat)
save(ENSAT, file = file.path(processed_path, "ENSAT_HLA_Clinical.RData"))

message("ENSAT Data Preparation Complete.")