#This script integrate individual H5 files into a single
#Seurat object. RStudio does not support parallelization
#with future, and this runs way faster through R CLI

rm(list = ls())
gc()

Sys.setenv(R_MAX_VSIZE = "100Gb")
options(future.globals.maxSize= 2000 * 1024^2)

library(tidyverse)
library(Seurat)
library(BiocGenerics)
library(patchwork)
library(hdf5r)
library(future)
plan("multisession", workers = 8)

setwd("./data/cohort2/ARDS-scRNAseq-unfiltered-samples/")
filenames <- dir()

samplename <- filenames[1] %>% str_extract("HS[:digit:]+")
ARDS.seu <-  Read10X_h5(filenames[1]) %>% CreateSeuratObject(project =samplename )

for (file in filenames[2:length(filenames)]) {
  samplename <- file %>% str_extract("HS[:digit:]+")
  tmp <- Read10X_h5(file) %>% CreateSeuratObject(project = samplename)
  ARDS.seu  <- merge(ARDS.seu, y = tmp)
  rm(tmp)
  gc()
}


ARDS.seu.list <- SplitObject(ARDS.seu,split.by = "orig.ident")
rm(ARDS.seu)
gc()

setwd("..") #Return to the cohort2 data folder

# normalize and identify variable features for each dataset independently
ARDS.seu.list <- lapply(X = ARDS.seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  Idents(x) <- "orig.ident"
  DefaultAssay(x) <- "RNA"
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  
  
  Idents(x) <- "major.label"
  DefaultAssay(x) <- "RNA"
  
  x <- subset(x, subset = 
                            nCount_RNA > 300 &
                            nCount_RNA < 30000 & 
                            percent.mt < 10)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ARDS.seu.list)

# integration -
ARDS.anchors <- FindIntegrationAnchors(object.list = ARDS.seu.list, anchor.features = features)
# save(ARDS.anchors, file = "ARDSintegrationanchors.Rda")
# load("ARDSintegrationanchors.Rda")

ARDS.combined <- IntegrateData(anchorset = ARDS.anchors)
#save(ARDS.combined, file = "cohort2.scSeq.Rda")
rm(ARDS.anchors)
rm(ARDS.seu.list)
gc()


#Annotate the imported cells using SingleR
library(SingleR)
library(BiocParallel)
MulticoreParam(6)

#This is a database of microarray experiments for cells exposted to a variety of stimuli
hpca.se <- celldex::HumanPrimaryCellAtlasData()

DefaultAssay(ARDS.combined) <- "RNA"
somecells <- GetAssayData(object = ARDS.combined, slot = "data")

cell.annot <- SingleR(test = somecells, ref = hpca.se,assay.type.test = 1,
                      labels = hpca.se$label.fine,BPPARAM = MulticoreParam(6))

ARDS.combined$singleR.label <- cell.annot$labels
ARDS.combined$major.label <- cell.annot$labels %>% str_remove(":.+")

rm(somecells)
rm(cell.annot)
rm(hpca.se)

save(ARDS.combined, file = "cohort2.scSeq.Rda")




