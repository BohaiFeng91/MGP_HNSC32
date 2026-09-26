library(MatchIt)
library(survival)
library(cobalt)

set.seed(123456)
dir.create("results/PSM", recursive = TRUE, showWarnings = FALSE)

# 1. Load TCGA data
# Required columns:
# ID, group, T_stage, N_stage, Resection_margin,
# Radiation_Therapy, OS_time, OS,
# PFI_time, PFI, DSS_time, DSS

dat <- read.csv("data/TCGA_PSM_input.csv")

dat$group <- factor(dat$group, levels = c("LR", "HR"))

# 2. Propensity score matching
# Use the same complete-case criteria as the original analysis
vars <- c(
  "group", "T_stage", "N_stage",
  "Resection_margin", "Radiation_Therapy"
)

dat <- dat[complete.cases(dat[, vars]), ]

psm <- matchit(
  group ~ T_stage + N_stage +
    Resection_margin + Radiation_Therapy,
  data = dat,
  method = "nearest",
  distance = "glm",
  estimand = "ATT",
  ratio = 1,
  replace = FALSE,
  caliper = 0.2,
  std.caliper = TRUE
)

# 3. Matching results and covariate balance
print(summary(psm, standardize = TRUE))

balance <- cobalt::bal.tab(psm, un = TRUE)

write.csv(
  as.data.frame(balance$Balance),
  "results/PSM/Covariate_balance.csv"
)

matched <- match.data(psm)

write.csv(
  matched,
  "results/PSM/Matched_patients.csv",
  row.names = FALSE
)

# 4. Survival analysis after matching
endpoints <- c("OS", "PFI", "DSS")

for (endpoint in endpoints) {
  
  time_col <- paste0(endpoint, "_time")
  
  sub <- matched[
    complete.cases(matched[, c(time_col, endpoint)]),
  ]
  
  formula <- as.formula(
    paste0("Surv(", time_col, ", ", endpoint, ") ~ group")
  )
  
  # Kaplan-Meier analysis
  km <- survfit(formula, data = sub)
  
  pdf(paste0("results/PSM/PSM_", endpoint, "_KM.pdf"))
  plot(
    km,
    col = c("blue", "red"),
    xlab = "Time (days)",
    ylab = "Survival probability",
    main = paste("PSM:", endpoint)
  )
  legend(
    "bottomleft",
    legend = levels(sub$group),
    col = c("blue", "red"),
    lty = 1
  )
  dev.off()
  
  # Log-rank test
  logrank <- survdiff(formula, data = sub)
  
  p_value <- pchisq(
    logrank$chisq,
    df = length(logrank$n) - 1,
    lower.tail = FALSE
  )
  
  # Cox model with matched-pair robust variance
  cox <- coxph(
    update(formula, . ~ . + cluster(subclass)),
    data = sub
  )
  
  result <- data.frame(
    Endpoint = endpoint,
    N = nrow(sub),
    HR = exp(coef(cox))[1],
    Lower95CI = exp(confint(cox))[1, 1],
    Upper95CI = exp(confint(cox))[1, 2],
    Cox_P = summary(cox)$coefficients[1, "Pr(>|z|)"],
    Logrank_P = p_value
  )
  
  write.csv(
    result,
    paste0("results/PSM/PSM_", endpoint, "_Cox.csv"),
    row.names = FALSE
  )
  
  print(result)
}

saveRDS(psm, "results/PSM/PSM_model.rds")