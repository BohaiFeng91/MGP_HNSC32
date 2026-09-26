library(survival)

set.seed(123456)
B <- 1000

# dat: PFI_time, PFI and the top 100 gene expression values
genes <- read.csv("results/RSF_top100_genes.csv")$Gene
dat <- read.csv("data/TCGA_model_data.csv")
dat <- na.omit(dat[, c("PFI_time", "PFI", genes)])

formula <- as.formula(
  paste("Surv(PFI_time, PFI) ~",
        paste(genes, collapse = " + "))
)

# Original model
fit <- coxph(formula, data = dat)
model <- step(fit, direction = "both", k = 2, trace = 0)
final_genes <- names(coef(model))

# Bootstrap selection
selected <- vector("list", B)

for (i in seq_len(B)) {
  boot <- dat[sample(nrow(dat), replace = TRUE), ]
  
  result <- try({
    fit <- coxph(formula, data = boot)
    names(coef(step(fit, direction = "both",
                    k = 2, trace = 0)))
  }, silent = TRUE)
  
  if (!inherits(result, "try-error"))
    selected[[i]] <- result
}

# Selection frequency among successful iterations
selected <- Filter(Negate(is.null), selected)
stopifnot(length(selected) > 0)

frequency <- sapply(
  genes,
  function(g) mean(vapply(selected, function(x) g %in% x,
                          logical(1)))
)

result <- data.frame(
  Gene = genes,
  Selection_frequency = frequency,
  Final_model = genes %in% final_genes
)

dir.create("results", showWarnings = FALSE)

write.csv(
  result[order(-result$Selection_frequency), ],
  "results/Bootstrap_selection_frequency.csv",
  row.names = FALSE
)