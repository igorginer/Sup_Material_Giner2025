# -------------------------------------------------------------------------
# Script: 05_Survival_Analysis.R
# Description: Kaplan-Meier analysis for Groups, Stages, and Genes.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(survival, survminer, tidyverse, here)

# Load Data
load(here("03_ENSAT_ACC", "processed_data", "ENSAT_HLA_Clinical.RData"))
meta <- ENSAT[[2]]
gexp <- as.data.frame(t(ENSAT[[1]]))

output_dir <- here("03_ENSAT_ACC", "results", "survival")
if(!dir.exists(output_dir)) dir.create(output_dir)

# 1. C1A vs C1B
fit_c1 <- survfit(Surv(OS.time, OS) ~ Consensus.clustering.K2, data = meta)
pdf(file.path(output_dir, "KM_C1A_C1B.pdf"))
print(ggsurvplot(fit_c1, data = meta, pval = TRUE, palette = c("#e41a1c", "#377eb8")))
dev.off()

# 2. Stage (1/2 vs 3/4)
meta$Stage_Group <- ifelse(meta$ENSAT.staging %in% c(1,2), "Early", "Late")
fit_stg <- survfit(Surv(OS.time, OS) ~ Stage_Group, data = meta)
pdf(file.path(output_dir, "KM_Stages.pdf"))
print(ggsurvplot(fit_stg, data = meta, pval = TRUE, palette = "jco"))
dev.off()

# 3. Gene Loop
genes <- colnames(gexp)
res_genes <- data.frame()

for(g in genes) {
  df <- meta
  df$Expr <- gexp[[g]][match(rownames(df), rownames(gexp))]
  med <- median(df$Expr, na.rm=T)
  df$Strata <- ifelse(df$Expr > med, "High", "Low")
  
  fit <- survfit(Surv(OS.time, OS) ~ Strata, data = df)
  pval <- surv_pvalue(fit)$pval
  
  res_genes <- rbind(res_genes, data.frame(Gene=g, P_Val=pval))
  
  # Save plot
  pdf(file.path(output_dir, paste0("KM_", g, ".pdf")))
  print(ggsurvplot(fit, data=df, title=g, pval=TRUE))
  dev.off()
}

write.csv(res_genes, file.path(output_dir, "Survival_LogRank_Genes.csv"))