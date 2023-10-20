# CD4T Expression Signature from Grady Dataset

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(glmnet)
library(matrixStats)

COMMON_GENES <- as.character(read.table("results/universe.txt")$V1)
FILTER_ADDL <- 0.0015 #additional filter for LASSO-selected features

## Discovery Dataset w/ Covariates
expr_grady <- load_grady_expr()
expr_grady <- expr_grady[rownames(expr_grady) %in% COMMON_GENES, ]

covars_grady <- load_grady_covar()
covars_grady <- covars_grady[ , c("SampleID","CD4T")]

covars_grady <- subset(covars_grady, SampleID %in% colnames(expr_grady))
expr_grady <- expr_grady[ , colnames(expr_grady) %in% covars_grady$SampleID]
expr_grady <- expr_grady[ , match(covars_grady$SampleID, colnames(expr_grady))]

## Signature Discovery:
stopifnot(identical(colnames(expr_grady), covars_grady$SampleID))
modLasso <- glmnet(
  t(expr_grady), 
  MinMaxScale(covars_grady$CD4T), 
  family = gaussian, 
  alpha = 1
)
plot(modLasso)
lambda_best <- min(modLasso$lambda)
lambda_best #0.0009712538

## Extract Coefficients:
coefs <- coef(modLasso, s=lambda_best)
coefs <- data.frame(
  Gene = coefs@Dimnames[[1]][1+coefs@i],
  Coef = coefs@x
)
# write.csv(coefs, "results/lasso_model_coefs.csv", row.names=FALSE, quote=FALSE)
