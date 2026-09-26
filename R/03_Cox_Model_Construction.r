library(survival)

# Load the top 100 RSF-ranked genes
genes <- read.csv("results/RSF_top100_genes.csv")$Gene

# Load H&E-predicted gene expression and PFI data
exp <- read.csv("data/scores_DL_TCGA.csv")
sur <- read.csv("data/TCGA_PFI.csv")

dat <- merge(sur, exp, by = "ID")
dat <- na.omit(dat[, c("PFI_time", "PFI", genes)])

# Full multivariable Cox model
full_model <- coxph(
  as.formula(paste(
    "Surv(PFI_time, PFI) ~",
    paste(genes, collapse = " + ")
  )),
  data = dat
)

# Bidirectional stepwise selection using AIC
step_model <- step(
  full_model,
  direction = "both",
  k = 2,
  trace = 0
)

# Export selected genes and coefficients
gene_coef <- data.frame(
  Gene = names(coef(step_model)),
  Coef = as.numeric(coef(step_model))
)

dir.create("models", showWarnings = FALSE)

write.csv(
  gene_coef,
  "models/MGP_HNSC32_coefficients.csv",
  row.names = FALSE
)

saveRDS(step_model, "models/MGP_HNSC32_Cox.rds")