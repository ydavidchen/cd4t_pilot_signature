# Identify Common Genes shared among datasets

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(matrixStats)

## Intersect Discovery & Application sets:
expr_grady <- load_grady_expr()

expr_hiv <- load_hiv_expr()

mean(rownames(expr_grady) %in% rownames(expr_hiv)) #0.8354304
mean(rownames(expr_hiv) %in% rownames(expr_grady)) #0.7366166

COMMON_GENES <- intersect(rownames(expr_grady), rownames(expr_hiv))
length(COMMON_GENES) #12549

# write.table(COMMON_GENES, "results/universe.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Gene correlation between datasets:
expr_grady <- expr_grady[COMMON_GENES, ]
expr_hiv <- expr_hiv[COMMON_GENES, ]

avg_expr <- merge(
  data.frame(Gene = rownames(expr_grady), Discovery = rowMeans2(expr_grady)),
  data.frame(Gene = rownames(expr_hiv), Application = rowMeans2(expr_hiv)),
  by = "Gene"
)

ggplot(avg_expr, aes(Discovery, Application)) +
  geom_point() + 
  geom_smooth(method="lm") +
  THEME_SCATTER

cor.test(avg_expr$Discovery, avg_expr$Application)
# 95% CI 0.2198556 0.2528928
# cor 0.2364425
