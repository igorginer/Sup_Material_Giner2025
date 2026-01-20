# -------------------------------------------------------------------------
# Script: 06_Clinical_Associations.R
# Description: Associations between Steroid Phenotype and Clinical/TME features.
# -------------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, ggpubr, here, openxlsx)

load(here("data", "coldataACC.RData"))
load(here("data", "tcgaACC_pre_processed.RData"))
# Assuming external purity file exists from script 15
# If not, use thorsson file or data_purity.RData provided in data folder
load(here("data", "data_purity.RData")) 

output_dir <- here("02_TCGA_ACC", "results", "clinical")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Merge Data
coldataACC$barcode <- substr(coldataACC$barcode, 1, 15)
data_purity$barcode <- substr(data_purity$barcode, 1, 15)

df <- coldataACC %>%
  left_join(data_purity, by="barcode") %>%
  filter(!is.na(steroid))

# -------------------------------------------------------------------------
# 1. Quantitative Tests (Wilcoxon) - Purity/Stromal
# -------------------------------------------------------------------------
vars_to_test <- c("paper_purity", "Stromal.Fraction") # Ensure column names match loaded data

res_quanti <- data.frame()

for(v in vars_to_test) {
  if(v %in% colnames(df)) {
    stat <- wilcox.test(as.formula(paste(v, "~ steroid")), data = df)
    
    p <- ggboxplot(df, x = "steroid", y = v, 
                   color = "steroid", palette = c("#e41a1c", "#377eb8"),
                   add = "jitter", title = paste(v, "by Steroid")) +
      stat_pvalue_manual(stat, label = "p = {p.signif}")
    
    ggsave(file.path(output_dir, paste0("Boxplot_", v, ".pdf")), p)
    
    res_quanti <- rbind(res_quanti, data.frame(Var = v, P_Val = stat$p.value))
  }
}
write.csv(res_quanti, file.path(output_dir, "Quantitative_Tests.csv"))

# -------------------------------------------------------------------------
# 2. Qualitative Tests (Fisher)
# -------------------------------------------------------------------------
vars_quali <- c("tumor_stage", "cortisol.excess", "immune.subtype")
res_quali <- data.frame()

for(v in vars_quali) {
  tbl <- table(df[[v]], df$steroid)
  ft <- fisher.test(tbl, simulate.p.value = TRUE) # Simulate for large tables
  
  res_quali <- rbind(res_quali, data.frame(Var = v, P_Val = ft$p.value))
}
write.csv(res_quali, file.path(output_dir, "Qualitative_Tests.csv"))