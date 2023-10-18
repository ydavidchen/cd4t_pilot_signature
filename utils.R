# Utility Module

library(data.table)
library(ggplot2)
library(readxl)

STRNA <- c("", "--", "NA")

## Paths **masked upon git push**:
GRADY_DIR <- "*************MASKED*************"
HIV_DIR <- "*************MASKED*************"

# ---------------------  Data Loading Methods --------------------- 
load_grady_covar <- function(prefix=GRADY_DIR, fname="*****MASKED*****.csv") {
  covars <- readxl::read_excel(paste0(prefix,fname), sheet="covars", na=STRNA, trim_ws=TRUE)
  covars <- as.data.frame(covars)
  covars$SampleID <- paste0("s_", covars$SampleID)
  
  covars$RaceEth[! covars$RaceEth %in% c("African American","Caucasian") | is.na(covars$RaceEth)] <- "Mixed/Other/Unk"
  covars$RaceWhite <- covars$RaceEth == "Caucasian"
  
  cellperc <- read_excel(paste0(prefix, fname), sheet="immune", na=STRNA, trim_ws=TRUE)

  covars <- merge(covars, cellperc[ , -1], by="Accession")
  return(covars)
}


load_grady_expr <- function(prefix=GRADY_DIR, prim_key="ID", zscore=TRUE) {
  #'@description Aggregate, Scale, Matrixify expression data
  mat <- fread(paste0(prefix, "*****MASKED*****.csv"), data.table=FALSE, na.strings=STRNA)
  colnames(mat)[1] <- prim_key
  
  eAnnot <-  fread(paste0(prefix, "*****MASKED*****.txt"), data.table=FALSE, na.strings=STRNA)
  eAnnot <- eAnnot[ , c(prim_key,"Symbol")]
  eAnnot <- eAnnot[complete.cases(eAnnot), ]
  
  mat <- merge(mat, eAnnot, by=prim_key)
  mat <- mat[ , -1]
  mat <- aggregate(. ~ Symbol, data=mat, FUN=mean, na.rm=TRUE)
  
  rownames(mat) <- as.character(unlist(mat[ , 1]))
  mat <- data.matrix(mat[ , -1])
  if(zscore) mat <- scale(mat)
  
  return(mat)
}

load_hiv_covars <- function(prefix=HIV_DIR) {
  covars <- read.csv(paste0(prefix, "*****MASKED*****.csv"), strip.white=TRUE, stringsAsFactors=FALSE)
  covars$RaceWhite <- covars$Race == "White"
  covars$PercChangeCD4 <- 100 * (covars$ChangeCD4CountPerMl / covars$CD4_count_baseline)
  return(covars)
}

load_hiv_expr <- function(prefix=HIV_DIR, prim_key="ID", zscore=TRUE) {
  mat <- fread(paste0(prefix,"*****MASKED*****.txt"), data.table=FALSE, na.strings=STRNA)
  colnames(mat)[1] <- prim_key
  
  eAnnot <- fread(paste0(prefix,"*****MASKED*****.txt"), data.table=FALSE, na.strings=STRNA)
  eAnnot <- eAnnot[ , c(prim_key,"Symbol")]
  
  mat <- merge(mat, eAnnot, by=prim_key)
  mat <- mat[ , -1]
  mat <- aggregate(. ~ Symbol, data=mat, FUN=mean, na.rm=TRUE)
  
  rownames(mat) <- as.character(unlist(mat[ , 1]))
  mat <- data.matrix(mat[ , -1])
  if(zscore) mat <- scale(mat)
  
  return(mat)
}

# --------------------- General Methods --------------------- 
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
