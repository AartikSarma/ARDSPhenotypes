# ARDS Phenotypes - 3 - Single-cell differential expression.R
#
# Tests for genes differentially expressed between HyperARDS and HypoARDS
# within each of the five most common cell types in the tracheal aspirate
# single-cell RNA-seq (GSE200848), the five reported in Supplementary Data E5,
# using MAST's hurdle model. The model and settings are those of
# scripts/preprint_2022/ARDS Phenotypes - 4 - scMAST.R, adapted from
# https://github.com/RGLab/MAST/issues/147:
#
#   ~ phenotype + cngeneson + (1 | orig.ident)
#
# where cngeneson is the scaled cellular detection rate (the number of genes
# detected in each cell) and orig.ident is the subject. Expression is
# log2(count + 1) of the raw counts; genes are kept if detected in at least
# 10% of the cells of that type. The discrete component is fit with glmer
# (nAGQ = 0), without empirical Bayes shrinkage, and phenotype is tested with a
# likelihood ratio test on the combined hurdle model.
#
# Run from the repository root after script 0. Writes one file per cell type to
# data/results/ajrccm_2025/scseq_mast/, with one row per tested gene: hurdle
# log2 fold change of HyperARDS over HypoARDS with its 95% confidence interval,
# likelihood ratio p value and Benjamini-Hochberg FDR. The same results, without
# the confidence intervals and listed alphabetically by gene, form
#   data/results/ajrccm_2025/Supplementary Data E5 - scSeq MAST results.xlsx
# MAST reports no fold change or p value for a gene whose model failed to fit.

library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(MAST)

clean_directory <- file.path(here::here(), "data", "geo", "clean")
results_directory <- file.path(here::here(), "data", "results", "ajrccm_2025", "scseq_mast")
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

# MAST fits genes in parallel through parallel::mclapply.
options(mc.cores = max(1, parallel::detectCores() - 2))

scseq_seurat <- readRDS(file.path(clean_directory, "scseq_seurat.rds"))

top_cell_types <- c("Neutrophil", "Monocyte", "Monocyte-derived mac", "T cell", "Alveolar mac")
stopifnot(setequal(top_cell_types,
                   names(head(sort(table(scseq_seurat$topcells), decreasing = TRUE), 5))))

run_mast <- function(cell_type) {
  cell_sce <- as.SingleCellExperiment(
    subset(scseq_seurat, cells = colnames(scseq_seurat)[scseq_seurat$topcells == cell_type]),
    assay = "RNA")
  logcounts(cell_sce) <- log2(counts(cell_sce) + 1)
  cell_sca <- SceToSingleCellAssay(cell_sce)

  cellular_detection_rate <- colSums(assay(cell_sca) > 0)
  colData(cell_sca)$cngeneson <- scale(cellular_detection_rate)
  # HypoARDS is the reference level, so the phenotype coefficient is the log2
  # fold change of HyperARDS over HypoARDS.
  colData(cell_sca)$phenotype <- relevel(factor(colData(cell_sca)$phenotype), "HypoARDS")

  cell_sca_filtered <- filterLowExpressedGenes(cell_sca, threshold = 0.1)

  zlm_fit <- zlm(~ phenotype + cngeneson + (1 | orig.ident),
                 sca = cell_sca_filtered,
                 exprs_value = "logcounts",
                 method = "glmer",
                 ebayes = FALSE,
                 silent = TRUE,
                 fitArgsD = list(nAGQ = 0),
                 strictConvergence = FALSE)

  phenotype_summary <- summary(zlm_fit, doLRT = "phenotypeHyperARDS")$datatable %>%
    as_tibble() %>%
    filter(contrast == "phenotypeHyperARDS")

  hurdle_results <-
    phenotype_summary %>%
    filter(component == "H") %>%
    dplyr::select(gene_symbol = primerid, p_value = `Pr(>Chisq)`) %>%
    inner_join(phenotype_summary %>%
                 filter(component == "logFC") %>%
                 dplyr::select(gene_symbol = primerid, log2_fold_change = coef,
                               ci_low = ci.lo, ci_high = ci.hi),
               by = "gene_symbol") %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value) %>%
    mutate(cell_type = cell_type, cells = ncol(cell_sca_filtered), .before = 1)

  write_csv(hurdle_results,
            file.path(results_directory, paste0(str_replace_all(cell_type, " ", "_"), ".csv")))
  message(cell_type, ": ", ncol(cell_sca_filtered), " cells, ", nrow(hurdle_results),
          " genes tested, ", sum(hurdle_results$fdr < 0.1, na.rm = TRUE), " with FDR < 0.1")
  hurdle_results
}

mast_results <- map(set_names(top_cell_types), run_mast)

# Supplementary Data E5 ----------------------------------------------------

# Sheet for each cell type, in the published order. Excel limits sheet names
# to 31 characters.
supplementary_sheets <- c(
  "T cell" = "S5A - T cells",
  "Monocyte-derived mac" = "S5B - Monocyte-derived macs",
  "Alveolar mac" = "S5C - Alveolar macrophages",
  "Monocyte" = "S5D - Monocytes",
  "Neutrophil" = "S5E - Neutrophils"
)

mast_results[names(supplementary_sheets)] %>%
  map(~ .x %>%
        arrange(gene_symbol) %>%
        transmute(HGNC = gene_symbol, `Log2 Fold Change` = log2_fold_change,
                  `P-value` = p_value, FDR = fdr)) %>%
  set_names(supplementary_sheets) %>%
  writexl::write_xlsx(file.path(dirname(results_directory),
                                "Supplementary Data E5 - scSeq MAST results.xlsx"))
