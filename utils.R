# Utility Module

library(data.table)
library(ggplot2)
library(readxl)

## Paths **masked upon git push**:
GRADY_DIR <- "*********** MASKED ***********"
HIV_DIR <- "*********** MASKED ***********"

STRNA <- c("", "--", "NA", "n/a")
NORM_VAL_BOUNDS <- c(-5, 5)

# ---------------------  Data Loading Methods --------------------- 
load_grady_covar <- function(prefix=GRADY_DIR, fname="GSE72680_series_matrix.xlsx") {
  covars <- readxl::read_excel(paste0(prefix,fname), sheet="covars", na=STRNA, trim_ws=TRUE)
  covars <- as.data.frame(covars)
  covars$SampleID <- paste0("s_", covars$SampleID)
  
  covars$RaceEth[! covars$RaceEth %in% c("African American","Caucasian") | is.na(covars$RaceEth)] <- "Mixed/Other/Unk"
  covars$RaceWhite <- covars$RaceEth == "Caucasian"
  
  cellperc <- read_excel(paste0(prefix, fname), sheet="immune", na=STRNA, trim_ws=TRUE)

  covars <- merge(covars, cellperc[ , -1], by="Accession")
  return(covars)
}


load_grady_expr <- function(prefix=GRADY_DIR, prim_key="ID", constraints=NORM_VAL_BOUNDS) {
  #'@description Aggregate, Scale, Matrixify expression data
  mat <- fread(paste0(prefix, "GSE58137_full_expression.csv"), data.table=FALSE, na.strings=STRNA)
  colnames(mat)[1] <- prim_key
  
  eAnnot <-  fread(paste0(prefix, "GPL6947-13512.txt"), data.table=FALSE, na.strings=STRNA)
  eAnnot <- eAnnot[ , c(prim_key,"Symbol")]
  eAnnot <- eAnnot[complete.cases(eAnnot), ]
  
  mat <- merge(eAnnot, mat, by=prim_key)
  mat <- mat[ , -1]
  mat <- aggregate(. ~ Symbol, data=mat, FUN=mean, na.rm=TRUE)
  
  rownames(mat) <- as.character(unlist(mat[ , 1]))
  
  mat <- data.matrix(mat[ , -1])
  
  mat <- scale(mat)
  
  if(! is.null(constraints)) {
    mat[mat < constraints[1]] <- constraints[1]
    mat[mat > constraints[2]] <- constraints[2]
  }
  
  return(mat)
}

load_hiv_covars <- function(prefix=HIV_DIR) {
  covars <- read.csv(paste0(prefix, "GSE19087_series_matrix.csv"), strip.white=TRUE, stringsAsFactors=FALSE)
  covars$RaceWhite <- covars$Race == "White"
  covars$PercChangeCD4 <- 100 * (covars$ChangeCD4CountPerMl / covars$CD4_count_baseline)
  return(covars)
}

load_hiv_expr <- function(prefix=HIV_DIR, prim_key="ID", constraints=NORM_VAL_BOUNDS) {
  mat <- fread(paste0(prefix,"GSE19087_processed_expr.txt"), data.table=FALSE, na.strings=STRNA)
  colnames(mat)[1] <- prim_key
  
  eAnnot <- fread(paste0(prefix,"GPL6884-11607.txt"), data.table=FALSE, na.strings=STRNA)
  eAnnot <- eAnnot[ , c(prim_key,"Symbol")]
  
  mat <- merge(mat, eAnnot, by=prim_key)
  mat <- mat[ , -1]
  mat <- aggregate(. ~ Symbol, data=mat, FUN=mean, na.rm=TRUE)
  
  rownames(mat) <- as.character(unlist(mat[ , 1]))
  
  mat <- data.matrix(mat[ , -1])
  
  mat <- scale(mat)
  
  if(! is.null(constraints)) {
    mat[mat < constraints[1]] <- constraints[1]
    mat[mat > constraints[2]] <- constraints[2]
  }
  
  return(mat)
}

# --------------------- General Methods --------------------- 
MinMaxScale <- function(x, na.rm=TRUE) return((x- min(x)) / (max(x)-min(x)))
undo_minmax <- function(x, ref_dat) return( (x - min(ref_dat)) / (max(ref_dat) - min(ref_dat)) )

eval_regres <- function(preds, truths) {
  #'@description Evaluates reression outputs
  rmse <- sqrt(mean( (preds - truths)^2 ))
  pcc <- cor(preds, truths)
  return( c(RMSE=rmse, Pearson=pcc, R2=pcc^2) )
}

custom_hier_clust <- function(tMat, num_cl=2, method_hc=CLUSTPARAM[1], method_dist=CLUSTPARAM[2]) {
  #'@param tMat Expression matrix where row=Samples(observations), col=genes(features)
  #'@references Chen 2022 Med.Res.Arch.
  myDist <- dist(tMat, method=method_dist)
  myHcl <- hclust(myDist, method=method_hc)
  
  plot(myHcl, labels=FALSE)
  rect.hclust(myHcl, k=num_cl)
  
  membership <- data.frame(cutree(myHcl, k=num_cl))
  colnames(membership) <- "Cluster"
  membership$Cluster <- paste0("Cluster", membership$Cluster)
  
  return(membership)
}

# --------------------- Plotting --------------------- 
BINARY_COLORS <- c(Yes="black", No="lightgray")

HEAT_COLS_EXPR <- colorRampPalette(c("blue","lightgray","red"))(1024)

THEME_SCATTER <- theme_bw() + 
  theme(axis.text.x=element_text(size=20,color="black"), axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
        panel.spacing = unit(0.5, "lines"), axis.line=element_line(color="black"),
        strip.text.x=element_text(size=20,color="black",face="bold"), strip.background=element_rect(fill="gray95"),
        legend.position="top", legend.text=element_text(size=15,color="black"))

THEME_BOX <- theme_bw() +
  theme(axis.text.x=element_text(size=10,color="black", angle=90), axis.text.y=element_text(size=21,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=20,color="black"),
        strip.text.x=element_text(size=20,color="black",face="bold"), strip.background=element_rect(fill="gray95"),
        panel.border = element_blank(), axis.line=element_line(color="black"),
        legend.position="top", legend.title=element_text(size=20), legend.text=element_text(size=15,color="black"))

THEME_BAR <- theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=15,color="black"),
        axis.title.x=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
        panel.border = element_blank(), axis.line=element_line(color="black"),
        legend.position="top", legend.title=element_text(size=20), legend.text=element_text(size=15,color="black"))
