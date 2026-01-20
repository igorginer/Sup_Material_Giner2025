# -------------------------------------------------------------------------
# Script: 04_Survival_Analysis.R
# Description: Kaplan-Meier and Cox Regression analyses.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(survival, survminer, tidyverse, matrixStats, here)

# Load Data
load(here("02_TCGA_ACC", "processed_data", "rnaseq_hla.RData"))
load(here("02_TCGA_ACC", "processed_data", "coldataACC_Clean.RData"))
survival_data <- read_tsv(here("data", "TCGA-ACC.survival.tsv.gz"))

output_dir <- here("02_TCGA_ACC", "results", "survival")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# Data Prep
survival_data$sample <- substr(survival_data$sample, 1, 15)
coldataACC$barcode <- substr(coldataACC$barcode, 1, 15)

# Merge
surv_df <- survival_data %>%
  inner_join(coldataACC, by = c("sample" = "barcode")) %>%
  filter(!is.na(OS.time), !is.na(OS))

# -------------------------------------------------------------------------
# 1. Clinical Groups KM (HSP vs LSP / C1A vs C1B)
# -------------------------------------------------------------------------
# HSP/LSP
fit_steroid <- survfit(Surv(OS.time, OS) ~ steroid, data = surv_df)
pdf(file.path(output_dir, "KM_Steroid.pdf"))
print(ggsurvplot(fit_steroid, data = surv_df, pval = TRUE, risk.table = TRUE, palette = c("#e41a1c", "#377eb8")))
dev.off()

# Stage Grouping (I/II vs III/IV)
surv_df <- surv_df %>%
  mutate(Stage_Group = case_when(
    tumor_stage %in% c("stage i", "stage ii") ~ "Early",
    tumor_stage %in% c("stage iii", "stage iv") ~ "Late",
    TRUE ~ NA_character_
  )) %>% filter(!is.na(Stage_Group))

fit_stage <- survfit(Surv(OS.time, OS) ~ Stage_Group, data = surv_df)
pdf(file.path(output_dir, "KM_Stage.pdf"))
print(ggsurvplot(fit_stage, data = surv_df, pval = TRUE, palette = "jco"))
dev.off()

# -------------------------------------------------------------------------
# 2. Gene-Level KM (Loop)
# -------------------------------------------------------------------------
# Prepare HLA Expression
colnames(rnaseq_hla)[-c(1,2)] <- substr(colnames(rnaseq_hla)[-c(1,2)], 1, 15)
expr_t <- rnaseq_hla %>% 
  select(-gene_id) %>%
  column_to_rownames("gene_name") %>% 
  t() %>% as.data.frame() %>%
  rownames_to_column("sample")

surv_gene_df <- surv_df %>% inner_join(expr_t, by = "sample")

res_gene_km <- data.frame()
genes <- colnames(expr_t)[-1]

for(g in genes) {
  med <- median(surv_gene_df[[g]], na.rm=T)
  surv_gene_df$Strata <- ifelse(surv_gene_df[[g]] > med, "High", "Low")
  
  fit <- survfit(Surv(OS.time, OS) ~ Strata, data = surv_gene_df)
  pval <- surv_pvalue(fit)$pval
  
  res_gene_km <- rbind(res_gene_km, data.frame(Gene = g, P_LogRank = pval))
}
write.csv(res_gene_km, file.path(output_dir, "HLA_Gene_LogRank_Results.csv"))

# -------------------------------------------------------------------------
# 3. Signature Score Analysis
# -------------------------------------------------------------------------
# Define Signature
sig_genes <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", 
               "HLA-DPA1", "HLA-DPB1", "HLA-DPB2", "HLA-DQA1", "HLA-DQA2",
               "HLA-DQB1", "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", 
               "HLA-DRB6", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB")

# Calculate Score
valid_genes <- intersect(sig_genes, colnames(surv_gene_df))
surv_gene_df$Signature_Score <- rowMeans(scale(surv_gene_df[, valid_genes]), na.rm=T)

# 4 Groups: Steroid + Score
med_sig <- median(surv_gene_df$Signature_Score, na.rm=T)
surv_gene_df$Sig_Status <- ifelse(surv_gene_df$Signature_Score > med_sig, "ScoreHigh", "ScoreLow")
surv_gene_df$Combined <- paste(surv_gene_df$steroid, surv_gene_df$Sig_Status, sep = "_")

fit_comb <- survfit(Surv(OS.time, OS) ~ Combined, data = surv_gene_df)
pdf(file.path(output_dir, "KM_4_Groups_Steroid_Signature.pdf"), width = 10)
print(ggsurvplot(fit_comb, data = surv_gene_df, pval = TRUE, risk.table=TRUE))
dev.off()

# -------------------------------------------------------------------------
# 4. Multivariate Cox Regression
# -------------------------------------------------------------------------
# Adjusting for Stage and Steroid
cox_res <- data.frame()

for(g in genes) {
  # Univariate
  f_uni <- as.formula(paste("Surv(OS.time, OS) ~", g))
  fit_uni <- coxph(f_uni, data = surv_gene_df)
  
  # Multivariate
  f_multi <- as.formula(paste("Surv(OS.time, OS) ~", g, "+ tumor_stage + steroid"))
  fit_multi <- coxph(f_multi, data = surv_gene_df)
  
  s_multi <- summary(fit_multi)
  
  cox_res <- rbind(cox_res, data.frame(
    Gene = g,
    HR_Multi = s_multi$conf.int[1, "exp(coef)"],
    P_Multi = s_multi$coefficients[1, "Pr(>|z|)"]
  ))
}
write.csv(cox_res, file.path(output_dir, "Cox_Multivariate_Results.csv"))