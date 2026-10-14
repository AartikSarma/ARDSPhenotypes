# ARDS Phenotypes - 6 - Ligand-receptor analysis.R
#
# Infers ligand-receptor communication between cell types in the tracheal
# aspirate single-cell RNA-seq (GSE200848) with CellChat, separately in
# HyperARDS and HypoARDS, and lists the interactions whose ligand and receptor
# are both more highly expressed in one phenotype than in the other
# (Supplementary Data E6).
#
# Steps, using CellChat 1.5.0 and its human database (CellChatDB v1):
#   1. Keep eight cell types: T cells, monocytes, dendritic cells, neutrophils,
#      alveolar macrophages, monocyte-derived macrophages, NK cells, B cells.
#   2. For each phenotype: find over-expressed genes and interactions, compute
#      communication probabilities (CellChat defaults: trimean, 100
#      permutations, seed 1), drop interactions in cell types with fewer than
#      10 cells, and aggregate the network.
#   3. Merge the two phenotypes, test genes for differential expression
#      between them within each cell type (genes detected in at least 10% of
#      cells, log fold change threshold 0.1, no p value threshold), and map the
#      results onto the interactions.
#   4. Keep interactions with communication p < 0.05 whose ligand and receptor
#      are both up in that phenotype: in HyperARDS, ligand log fold change
#      >= 0.2 and receptor >= 0.1; in HypoARDS, both <= -0.1.
# Expression is the log-normalized "data" layer of the GEO object.
#
# CellChat 1.5.0 is installed from its GitHub tag (sqjin/CellChat@v1.5.0).
# The tag includes compiled x86_64 objects in src/; on Apple silicon, delete
# them and build the package from source.
#
# Run from the repository root after script 0. Writes to
# data/results/ajrccm_2025/:
#   Supplementary Data E6 - Ligand-receptor pairs.xlsx
#   cellchat_hyperards_hypoards.rds   merged CellChat object (HypoARDS, HyperARDS)

library(tidyverse)
library(Seurat)
library(CellChat)

clean_directory <- file.path(here::here(), "data", "geo", "clean")
results_directory <- file.path(here::here(), "data", "results", "ajrccm_2025")
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

stopifnot(packageVersion("CellChat") == "1.5.0")

cell_types <- c("T cell", "Monocyte", "Dendritic cell", "Neutrophil", "Alveolar mac",
                "Monocyte-derived mac", "NK cell", "B cell")
phenotypes <- c("HypoARDS", "HyperARDS")

future::plan("multisession", workers = 6)

scseq_seurat <- readRDS(file.path(clean_directory, "scseq_seurat.rds"))
scseq_seurat <- subset(scseq_seurat,
                       cells = colnames(scseq_seurat)[scseq_seurat$topcells %in% cell_types])

infer_communication <- function(phenotype) {
  phenotype_seurat <- subset(scseq_seurat,
                             cells = colnames(scseq_seurat)[scseq_seurat$phenotype == phenotype])
  # CellChat 1.5.0 reads a Seurat object with an accessor that SeuratObject 5
  # no longer supports, so it is given the expression matrix and cell metadata
  # that it would extract.
  cellchat <- createCellChat(LayerData(phenotype_seurat, assay = "RNA", layer = "data"),
                             meta = phenotype_seurat[[]], group.by = "topcells")
  cellchat@DB <- CellChatDB.human
  cellchat %>%
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb() %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet()
}

phenotype_cellchat <- map(set_names(phenotypes), infer_communication)
merged_cellchat <- mergeCellChat(phenotype_cellchat, add.names = phenotypes)

merged_cellchat <- identifyOverExpressedGenes(merged_cellchat, group.dataset = "datasets",
                                              pos.dataset = "HyperARDS",
                                              features.name = "HyperARDS", only.pos = FALSE,
                                              thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
mapped_interactions <- netMappingDEG(merged_cellchat, features.name = "HyperARDS")

phenotype_interactions <- list(
  "6A - Hyperinflammatory ARDS" =
    subsetCommunication(merged_cellchat, net = mapped_interactions, datasets = "HyperARDS",
                        ligand.logFC = 0.2, receptor.logFC = 0.1),
  "6B - Hypoinflammatory ARDS" =
    subsetCommunication(merged_cellchat, net = mapped_interactions, datasets = "HypoARDS",
                        ligand.logFC = -0.1, receptor.logFC = -0.1)
) %>%
  map(as_tibble)

writexl::write_xlsx(phenotype_interactions,
                    file.path(results_directory, "Supplementary Data E6 - Ligand-receptor pairs.xlsx"))

# The merged object is saved for the figures (script 7), which need the
# interaction networks of both phenotypes.
saveRDS(merged_cellchat, file.path(results_directory, "cellchat_hyperards_hypoards.rds"))
