#rm(list = ls())
setwd("~/Research/ARDSPhenotypes/")
gc()
Sys.setenv(R_MAX_VSIZE = "100Gb")

library(tidyverse)
library(Seurat)
library(BiocGenerics)
library(patchwork)
library(hdf5r)
library(future)
plan("multicore", workers = 8)

load("./data/cohort2/cohort2.scSeq.Rda")
DefaultAssay(ARDS.combined) <- "integrated"

ARDS.combined <- subset(ARDS.combined, orig.ident %in% paste("HS", c(16, 21,347,358,383,408, 447, 449), sep = ""))

#D0 biomarker assignment - 0.5
HyperARDS.samples <- paste("HS", c(16, 408), sep = "")
HypoARDS.samples <- paste("HS", c(347, 358, 21, 383, 447, 449), sep = "")

Idents(ARDS.combined) <- "orig.ident"
DefaultAssay(ARDS.combined) <- "RNA"
ARDS.combined[["percent.mt"]] <- PercentageFeatureSet(ARDS.combined, pattern = "^MT-")


Idents(ARDS.combined) <- "major.label"
DefaultAssay(ARDS.combined) <- "RNA"

ARDS.combined <- subset(ARDS.combined, subset = 
                          nCount_RNA > 300 &
                          nCount_RNA < 30000 & 
                          percent.mt < 10)

ARDS.combined$phenotype <- 
  case_when(
    ARDS.combined$orig.ident %in% HyperARDS.samples ~ "HyperARDS", 
    ARDS.combined$orig.ident %in% HypoARDS.samples ~ "HypoARDS", 
  )

ARDS.combined <- subset(ARDS.combined, phenotype %in% c("HyperARDS", "HypoARDS"))

phenotable <- cbind(ARDS.combined$orig.ident, ARDS.combined$phenotype) %>% 
  as.data.frame %>%
  distinct() %>% 
  dplyr::rename(
    subject = V1, 
    phenotype = V2
  ) %>%
  group_by(phenotype) %>%
  mutate(phenonum = row_number(), 
         phenonum = paste(phenotype, phenonum, sep = "")) %>%
  ungroup %>%
  dplyr::select(-phenotype)

ARDS.combined$subject <- ARDS.combined$orig.ident %>%
  enframe %>%
  dplyr::select(-name)  %>%
  dplyr::rename(subject = value) %>%
  left_join(phenotable) %>%
  dplyr::select(-subject) %>%
  unlist


# Run the standard workflow for visualization and clustering
DefaultAssay(ARDS.combined) <- "integrated"
ARDS.combined <- ScaleData(ARDS.combined, verbose = TRUE,vars.to.regress = "subject")
ARDS.combined <- RunPCA(ARDS.combined, npcs = 30 , verbose = FALSE)
ARDS.combined <- RunUMAP(ARDS.combined, reduction = "pca", dims = 1:30)
ARDS.combined <- FindNeighbors(ARDS.combined, reduction = "pca", dims = 1:30)
ARDS.combined <- FindClusters(ARDS.combined, resolution = 1)


DefaultAssay(ARDS.combined) <- "RNA"

Idents(ARDS.combined) <- "subject"
subj.umap <- DimPlot(ARDS.combined, reduction = "umap")
Idents(ARDS.combined) <- "phenotype"
pheno.umap <- DimPlot(ARDS.combined, reduction = "umap")
Idents(ARDS.combined) <- "seurat_clusters"
cluster.umap <- DimPlot(ARDS.combined, reduction = "umap", label = TRUE, repel = TRUE)

topcelltypes <- table(ARDS.combined$phenotype, ARDS.combined$singleR.label) %>%
  as.data.frame() %>%
  dplyr::rename(tx = Var1, celltype = Var2) %>%
  group_by(celltype) %>%
  filter(sum(Freq) > 100) %>%
  group_by(tx) %>%
  mutate(cellpct = Freq/sum(Freq) * 100) %>%
  ungroup %>%
  dplyr::select(celltype) %>%
  unique %>% 
  unlist

ARDS.combined$topcells <- case_when(ARDS.combined$major.label == "Neutrophil" ~ "Neutrophil", 
                                    ARDS.combined$major.label == "B_cell" ~ "B cell", 
                                    ARDS.combined$major.label == "DC" ~ "Dendritic cell", 
                                    ARDS.combined$major.label == "Epithelial_cells" ~ "Epithelial cell", 
                                    ARDS.combined$major.label == "HSC_-G-CSF" ~ "Hematopoietic stem cell",
                                    str_detect(ARDS.combined$singleR.label, "Macrophage:Alveolar") ~ "Alveolar mac", 
                                    str_detect(ARDS.combined$singleR.label, "Macrophage:monocyte-derived") ~ "Monocyte-derived mac", 
                                    ARDS.combined$major.label == "Monocyte" ~ "Monocyte",
                                    ARDS.combined$major.label == "NK_cell" ~ "NK cell",
                                    ARDS.combined$major.label == "T_cell" ~ "T cell", 
                                    TRUE ~ "Other")

Idents(ARDS.combined) <- "topcells"

save(ARDS.combined, file = "./data/cohort2/cohort2.scSeq.Rda")