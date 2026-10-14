#Use MAST to calculate differential expression for the most common cells in Seurat
#object. Future is not supported in RStudio - strongly recommend running this script
# from the R CLI. 

setwd("~/Research/ARDSPhenotypes/")
load("./data/cohort2/cohort2.scSeq.Rda")

setwd("output")
dir.create("scMAST")
setwd("scMAST")

library(tidyverse)
library(Seurat)
library(future)
library(magrittr)
library(SingleCellExperiment)
library(MAST)

Idents(ARDS.combined) <- "topcells"
plan(multicore, workers = 8)

ARDS.combined$topcells %>%
  table %>%
  sort(decreasing = TRUE) %>% 
  head(10) %>% names -> topcells

#MAST code adapted from: https://github.com/RGLab/MAST/issues/147

for(cell in topcells){
  
  cell.sce <- as.SingleCellExperiment(subset(ARDS.combined, topcells == cell), assay = "RNA")
  
  logcounts(cell.sce) <- log2(counts(cell.sce) + 1) #log2 transform the raw data**
  cell.sca <- SceToSingleCellAssay(cell.sce) 
  
  #Calculate the cellular detection rate
  cdr <-colSums(assay(cell.sca)>0)
  colData(cell.sca)$cngeneson <- scale(cdr)
  
  #Set the reference condition to HypoARDS
  cond <-factor(colData(cell.sca)$phenotype)
  cond <-relevel(cond,"HypoARDS")
  colData(cell.sca)$phenotype <-cond
  
  #Filter genes to those detected in at least 10% of cells
  cell.sca_filtered<-filterLowExpressedGenes(cell.sca,threshold=0.1)
  
  #Fit the mixed effects model:
  zlm_test <- zlm(~phenotype + # Fixed effect for phenotype
                    cngeneson + # Fixed effect for CDR
                    (1 | orig.ident), #Random effect for subject
                  sca = cell.sca_filtered,
                  exprs_value = 'logcounts',
                  method="glmer", 
                  ebayes=FALSE,
                  silent=T, 
                  fitArgsD = list(nAGQ = 0),
                  strictConvergence = FALSE)
  
  #Likelihood ratio test
  summaryCondtest <- summary(zlm_test, doLRT='phenotypeHyperARDS') 
  summaryDttest <- summaryCondtest$datatable
  
  fcHurdletest <- merge(summaryDttest[contrast=='phenotypeHyperARDS' & component=='H',.(primerid, `Pr(>Chisq)`)], 
                        summaryDttest[contrast=='phenotypeHyperARDS' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], 
                        by='primerid')
  
  fcHurdletest[,fdr:=p.adjust(`Pr(>Chisq)`,method ="BH")]
  
  saveRDS(fcHurdletest, file = paste(cell, "-5pct MAST results.Rds"))

}

