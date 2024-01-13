# Signature Validation by K-Fold
# Methods & Supplemental Results

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(glmnet)

COMMON_GENES <- as.character(read.table("results/universe.txt")$V1)
LAMBDA_BEST <- 0.0009712538 #best model hyperparam

wrapper_kfold <- function(data, label, nfolds=10) {
  #'@description Wrapper to run K-fold for evaluation
  #'@param data Input Features where Column=features (genes), Row=samples
  #'@param label Regression label
  #'@param n_folds Number of folds

  folds <- caret::createFolds(label, nfolds)
  
  res <- NULL
  for (k in 1:nfolds) {
    cat("Iteration: ", k, " of ", nfolds, "\n")
    
    ## Partition data into training and hold-out sets:
    idx_train <- unlist(folds[-k])
    idx_test <- folds[[k]]
    
    ## Modeling:
    model <- glmnet(data[idx_train, ], label[idx_train], family=gaussian, alpha=1, lambda=LAMBDA_BEST)
    pred_test <- predict(model, newx=data[idx_test, ])
    
    ## Evaluation:
    res <- rbind(res, eval_regres(pred_test, label[idx_test]))
  }
  rownames(res) <- paste0("Fold_", 1:nfolds)
  
  return(as.data.frame(res))
}

## Load data (same as before):
expr_grady <- load_grady_expr()
expr_grady <- expr_grady[rownames(expr_grady) %in% COMMON_GENES, ]

covars_grady <- load_grady_covar()
covars_grady <- covars_grady[ , c("SampleID","CD4T")]

covars_grady <- subset(covars_grady, SampleID %in% colnames(expr_grady))
expr_grady <- expr_grady[ , colnames(expr_grady) %in% covars_grady$SampleID]
expr_grady <- expr_grady[ , match(covars_grady$SampleID, colnames(expr_grady))]

## Implement K-fold for Evaluation: 
stopifnot( identical(covars_grady$SampleID, colnames(expr_grady)) ) #required checkpoint

cvRes <- wrapper_kfold(
  data = t(expr_grady), 
  label = MinMaxScale(covars_grady$CD4T), 
  nfolds = 10
)
cvRes$rawRMSE <- sapply(cvRes$RMSE, FUN=undo_minmax, ref_dat=covars_grady$CD4T) #empirically undo label transformation
cvRes <- rbind(cvRes, Avg=colMeans(cvRes))
cvRes

# write.csv(cvRes, "results/kfold_val.csv", row.names=TRUE, quote=FALSE)

