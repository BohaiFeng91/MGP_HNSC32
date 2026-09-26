library(survival)
library(randomSurvivalForest)

# RSF model
set.seed(123456)
fit <- rsf(
  Surv(PFI_time, PFI) ~ .,
  data = mydata,
  ntree = 1000,
  seed = 123456
)

# Gene importance
imp <- fit$importance[1, ]
imp <- sort(imp, decreasing = TRUE)

# Top 100 genes
top100 <- head(imp, 100)

write.csv(
  data.frame(Gene = names(top100),
             Importance = as.numeric(top100)),
  "RSF_top100_genes.csv",
  row.names = FALSE
)