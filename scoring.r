library(randomForestSRC)
score_with_models <- function(DL_new, rf_models) {
  ids <- DL_new$ID
  rownames(DL_new) <- ids
  X_all <- DL_new[, -which(colnames(DL_new) == "ID"), drop = FALSE]
  
  scores <- data.frame(ID = ids, stringsAsFactors = FALSE)
  
  for (gene_target in names(rf_models)) {
    info <- rf_models[[gene_target]]
    model <- info$model
    feats <- info$features
    
    # 准备新数据特征矩阵
    missing_feats <- setdiff(feats, colnames(X_all))
    X_new <- X_all
    if (length(missing_feats) > 0) {
      for (mf in missing_feats) X_new[, mf] <- NA
    }
    X_new_sel <- as.data.frame(X_new[, feats, drop = FALSE], check.names = FALSE)
    
    # 预测
    pred_out <- tryCatch({
      predict(model, newdata = X_new_sel)
    }, error = function(e) NULL)
    
    preds <- rep(NA_real_, nrow(X_new_sel))
    if (!is.null(pred_out)) preds <- as.numeric(pred_out$predicted)
    
    scores[[gene_target]] <- preds
  }
  
  # Z-score标准化
  gene_cols <- setdiff(colnames(scores), "ID")
  for (g in gene_cols) {
    scores[[g]] <- as.numeric(scale(scores[[g]], center = TRUE, scale = TRUE))
  }
  
  return(scores)
}
rf_models <- readRDS("rf_models.rds")


inputPath <- "features.csv"
DL<-read.csv(inputPath)
scores_DL <- score_with_models(DL, rf_models)
write.csv(scores_DL,"scores.csv",row.names = FALSE)
