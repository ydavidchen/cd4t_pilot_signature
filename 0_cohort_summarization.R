# Summarization of Cohort Information & Other
# Included in revised manuscript as new Supplemental Table 1

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(tableone)

covars_grady <- load_grady_covar()
covars_grady$Cohort <- "Discovery"

expr_grady <- load_grady_expr()
expr_grady <- expr_grady[ , colnames(expr_grady) %in% covars_grady$SampleID]

covars_grady <- subset(covars_grady, SampleID %in% colnames(expr_grady))
rm(expr_grady)

covars_hiv <- load_hiv_covars()
covars_hiv$Sex <- "Male"
covars_hiv$Cohort <- "Application"

## TODO: Stack 2 cohorts
COLS <- c("Accession","Cohort", "Age","Sex","RaceWhite")
meta <- rbind(covars_grady[,COLS], covars_hiv[,COLS])

t1 <- CreateTableOne(
  COLS[c(3,4,5)],
  strata = "Cohort",
  data = meta, 
  test = FALSE
)

y <- print(t1, showAllLevels=TRUE)
write.csv(y, "results/tableone.csv")
