library(randomForestSRC)
library(dplyr)

# 1. Load TCGA expression and pathology features
exp <- read.table("data/mRNA_FPKM.txt",
                 sep = "\t", header = TRUE,row.names = 1,
                 check.names = FALSE)
storage.mode(exp) <- "numeric"

features <- read.csv("data/features.csv")
features <- distinct(features, ID, .keep_all = TRUE)
rownames(features) <- features$ID

ids <- intersect(colnames(exp), features$ID)
X <- as.matrix(features[ids, -which(names(features) == "ID")])
storage.mode(X) <- "numeric"
exp <- exp[, ids, drop = FALSE]

# 2. Patient-level training/validation split
# Each ID must correspond to one patient.
set.seed(1000)
train_ids <- sample(ids, floor(length(ids) * 0.75))
test_ids <- setdiff(ids, train_ids)

X_train <- X[train_ids, , drop = FALSE]
X_test <- X[test_ids, , drop = FALSE]

# 3. Spearman screening and random forest regression
results <- list()

for (gene in rownames(exp)) {
  
  y_train <- as.numeric(exp[gene, train_ids])
  y_test <- as.numeric(exp[gene, test_ids])
  
  # Select pathology features in the training cohort
  pvals <- sapply(seq_len(ncol(X_train)), function(j) {
    suppressWarnings(
      tryCatch(
        cor.test(X_train[, j], y_train,
                 method = "spearman")$p.value,
        error = function(e) NA_real_
      )
    )
  })
  
  selected <- colnames(X_train)[
    !is.na(pvals) & pvals < 0.05
  ]
  if (length(selected) == 0) next
  
  # Gene-specific random forest
  set.seed(123456)
  model <- rfsrc(
    y ~ .,
    data = data.frame(
      y = y_train,
      X_train[, selected, drop = FALSE],
      check.names = FALSE
    ),
    ntree = 1000
  )
  
  # Validation
  pred <- predict(
    model,
    newdata = data.frame(
      X_test[, selected, drop = FALSE],
      check.names = FALSE
    )
  )$predicted
  
  test <- suppressWarnings(
    tryCatch(
      cor.test(pred, y_test, method = "spearman"),
      error = function(e) NULL
    )
  )
  if (is.null(test)) next
  
  results[[gene]] <- data.frame(
    Gene = gene,
    SpearmanRho = unname(test$estimate),
    p.value = test$p.value
  )
}

# 4. BH-FDR correction and final gene selection
results <- do.call(rbind, results)
results$FDR <- p.adjust(results$p.value, method = "BH")

genes <- results %>%
  filter(SpearmanRho >= 0.30, FDR <= 0.05)

dir.create("results", showWarnings = FALSE)

write.csv(results, "results/all_genes.csv", row.names = FALSE)
write.csv(genes, "results/morphology_genes.csv",
          row.names = FALSE)