library(survival)
library(maxstat)

# Load the final Cox model
model <- readRDS("models/MGP_HNSC32_Cox.rds")

# Load TCGA predicted gene expression
dat <- read.csv("data/scores_DL_TCGA.csv")

# Calculate relative risk
dat$riskScore <- as.numeric(
  predict(model, newdata = dat, type = "risk")
)

write.csv(
  dat[, c("ID", "riskScore")],
  "results/TCGA_risk_scores.csv",
  row.names = FALSE
)

# Load TCGA risk scores and survival data
score <- read.csv("results/TCGA_risk_scores.csv")
sur <- read.csv("data/TCGA_PFI.csv")

dat <- merge(score, sur, by = "ID")

# Determine the TCGA-derived cutoff
fit <- maxstat.test(
  Surv(PFI_time, PFI) ~ riskScore,
  data = dat,
  smethod = "LogRank",
  minprop = 0.2,
  maxprop = 0.8,
  pmethod = "exactGauss"
)

cutoff <- as.numeric(fit$estimate)

# Assign risk groups
dat$group <- ifelse(
  dat$riskScore > cutoff,
  "High", "Low"
)

# Save cutoff and patient risk groups
write.csv(
  data.frame(cutoff = cutoff),
  "models/MGP_HNSC32_cutoff.csv",
  row.names = FALSE
)

write.csv(
  dat,
  "results/TCGA_risk_groups.csv",
  row.names = FALSE
)

# Survival analysis
cox <- coxph(
  Surv(PFI_time, PFI) ~ group,
  data = transform(
    dat,
    group = relevel(factor(group), ref = "Low")
  )
)

summary(cox)

library(survival)

# Load the trained model and TCGA cutoff
model <- readRDS("models/MGP_HNSC32_Cox.rds")

cutoff <- read.csv(
  "models/MGP_HNSC32_cutoff.csv"
)$cutoff[1]

# Load HANCOCK predicted gene expression
exp <- read.csv("data/scores_DL_HANCOCK.csv")

# Predict risk scores without retraining
exp$riskScore <- as.numeric(
  predict(model, newdata = exp, type = "risk")
)

# Apply the unchanged TCGA cutoff
exp$group <- ifelse(
  exp$riskScore > cutoff,
  "High", "Low"
)

# Merge HANCOCK clinical data
clin <- read.csv("data/HANCOCK_clinical.csv")
dat <- merge(exp, clin, by = "ID")

dat$group <- factor(
  dat$group,
  levels = c("Low", "High")
)

# Kaplan-Meier analysis
km <- survfit(
  Surv(PFS_time, PFS) ~ group,
  data = dat
)

plot(km, xlab = "Days", ylab = "PFS probability")

# Grouped Cox regression
cox_group <- coxph(
  Surv(PFS_time, PFS) ~ group,
  data = dat
)

# Continuous risk score Cox regression
cox_continuous <- coxph(
  Surv(PFS_time, PFS) ~ riskScore,
  data = dat
)

summary(cox_group)
summary(cox_continuous)

# Export external validation results
write.csv(
  dat[, c("ID", "riskScore", "group",
          "PFS_time", "PFS")],
  "results/HANCOCK_validation.csv",
  row.names = FALSE
)