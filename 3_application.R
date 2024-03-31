# Signature Application via Heatmap Clustering

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)
library(tableone)

CLUSTPARAM <- c("ward.D2", "euclidean")
COEF_CUT <- 0.00125
CD4_CUT_DISC <- 0.25 #fraction
CD4_CUT_APPL <- 50 #percent

## LASSO Coefficients & Additinoal Filtering
coefs <- read.csv("results/lasso_model_coefs.csv", na.strings=STRNA)[-1 , ]
coefs$Keep <- ifelse(abs(coefs$Coef) > COEF_CUT, "Yes", "No")

ggplot(aes(reorder(Gene, Coef), Coef, fill=Keep), data=coefs[-1, ]) +
  geom_bar(stat="identity") +
  labs(x="Gene", y="LASSO Coefficient") +
  THEME_BAR

## Subset for downstream analyses & GOBP:
coefs <- subset(coefs, Keep=="Yes")
coefs$Keep <- NULL
coefs$Direction <- ifelse(coefs$Coef > 0, "Positive", "Negative")
table(coefs$Direction)
# write.csv(coefs, "results/lasso_filtered.csv", row.names = FALSE, quote=FALSE) #supp table online

# -------------------------- Visualization on Discovery Set --------------------------
## Grady Dataset - no need mutual subsetting for visualization
expr_grady <- load_grady_expr()
covars_grady <- load_grady_covar()

## Heatmap:
hm_grady_annot <- data.frame(
  row.names = covars_grady$SampleID, #or SampleID
  HighCD4T = ifelse(covars_grady$CD4T > CD4_CUT_DISC, "Yes", "No"),
  AgeOver35 = ifelse(covars_grady$Age > 35, "Yes", "No"),
  Male = ifelse(covars_grady$Sex=="Male", "Yes", "No"),
  Caucasian = ifelse(covars_grady$RaceWhite, "Yes", "No")
)

hm_grady_colors <- list()
hm_grady_colors[["HighCD4T"]] <- hm_grady_colors[["AgeOver35"]] <- hm_grady_colors[["Male"]] <- hm_grady_colors[["Caucasian"]] <- BINARY_COLORS

pheatmap(
  expr_grady[rownames(expr_grady) %in% coefs$Gene, ],
  show_rownames = FALSE,
  show_colnames = FALSE,
  cutree_cols = 2,
  clustering_method = CLUSTPARAM[1],
  annotation_colors = hm_grady_colors,
  annotation_col = hm_grady_annot,
  color = HEAT_COLS_EXPR,
  border_color = NA,
  fontsize = 12
)

res_cl_grady <- custom_hier_clust(t(expr_grady[rownames(expr_grady) %in% coefs$Gene, ]))
res_cl_grady$SampleID <- rownames(res_cl_grady)
res_cl_grady <- merge(res_cl_grady, covars_grady, by="SampleID")

ctab_grady <- table(
  Cluster2 = res_cl_grady$Cluster == "Cluster2",
  HighCD4T = res_cl_grady$CD4T > CD4_CUT_DISC
)
ctab_grady <- ctab_grady[c(2,1), c(2,1)]
ctab_grady

## Summarize by cluster - reviewer request:
t1_grady <- CreateTableOne(
  vars = c("Age","Sex","BMI","RaceWhite"),
  strata = "Cluster",
  data = res_cl_grady
)
print(t1_grady, showAllLevels=TRUE)

# -------------------------- Application Set --------------------------
expr_hiv <- load_hiv_expr()
covars_hiv <- load_hiv_covars()

hm_hiv_annot <- data.frame(
  row.names = covars_hiv$Accession,
  SAR = ifelse(covars_hiv$PercChangeCD4>=CD4_CUT_APPL, "Yes", "No"),
  Caucasian = ifelse(covars_hiv$RaceWhite, "Yes", "No"),
  AgeOver35 = ifelse(covars_hiv$Age >= 35, "Yes", "No") #median
)

hm_hiv_colors <- list()
hm_hiv_colors[["SAR"]] <- hm_hiv_colors[["AgeOver35"]] <- hm_hiv_colors[["Caucasian"]] <- BINARY_COLORS

pheatmap(
  expr_hiv[rownames(expr_hiv) %in% coefs$Gene, ],
  show_rownames = FALSE,
  show_colnames = FALSE,
  cutree_cols = 2,
  clustering_method = CLUSTPARAM[1],
  clustering_distance_rows = CLUSTPARAM[2],
  clustering_distance_cols = CLUSTPARAM[2],
  annotation_colors = hm_hiv_colors,
  annotation_col = hm_hiv_annot,
  color = HEAT_COLS_EXPR,
  border_color = NA,
  fontsize = 14
)

## Extract cluster membership for association/enrichment analysis:
res_cl <- custom_hier_clust(t(expr_hiv[rownames(expr_hiv) %in% coefs$Gene, ]))
res_cl$Accession <- rownames(res_cl)
res_cl <- merge(res_cl, covars_hiv, by="Accession")

## Fisher's test:
cTabUnivar <- table(
  Cluster1 = res_cl$Cluster=="Cluster1", 
  OverThresh = res_cl$PercChangeCD4 >= CD4_CUT_APPL
)
cTabUnivar <- cTabUnivar[c(2,1), c(2,1)]
cTabUnivar
fisher.test(cTabUnivar) 

## Summarize by cluster - reviewer request:
t1_h1v <- CreateTableOne(
  c("Age","RaceWhite", "CD4_count_baseline", "CD4_count_week48", "PercChangeCD4"),
  strata = "Cluster",
  data = res_cl,
  test = TRUE
)

print(t1_h1v, showAllLevels=TRUE)

## With covariate adjustment:
mlr <- glm(
  factor(Cluster=="Cluster1") ~ factor(PercChangeCD4 >= CD4_CUT_APPL) + Age + RaceWhite,
  data = res_cl,
  family = binomial
)
summary(mlr)
cbind(
  Estimate = exp(mlr$coefficients),
  exp(confint(mlr))
)

## Compare mean difference in %change in CD4T:
t.test(PercChangeCD4 ~ Cluster, data=res_cl)

## With covariate adjustment:
linreg <- glm(PercChangeCD4 ~ factor(Cluster=="Cluster1") + Age + RaceWhite, data=res_cl)
summary(linreg)
cbind(
  Estimate = linreg$coefficients,
  confint(linreg)
)
