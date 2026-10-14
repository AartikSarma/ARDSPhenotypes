# ARDS Phenotypes - 1 - Non-COVID bulk differential expression.R
#
# Tests for genes differentially expressed between HyperARDS and HypoARDS in
# tracheal aspirate bulk RNA-seq from the non-COVID cohort (GSE200847), one
# sample per subject, with limma-voom:
#
#   ~ rna_batch + phenotype + sex + age
#
# The steps follow rnaseq_to_limma() in tidyseq v0.1.2
# (https://github.com/AartikSarma/tidyseq/tree/v0.1.2), which fit the published
# model: genes with zero counts in every sample are removed, library sizes are
# TMM-normalized, genes are kept by filterByExpr() with an intercept-only
# design, and the model is fit with voomLmFit() followed by robust empirical
# Bayes moderation. Age enters unscaled, as in the published model.
#
# Age is top-coded at 90 (script 0), whereas the published model used exact
# ages, so results may differ subtly from those published. Additional data are
# available on reasonable request.
#
# Run from the repository root after script 0. Writes
#   data/results/ajrccm_2025/noncovid_bulk_hyper_vs_hypo.csv
# with one row per tested gene: log2 fold change of HyperARDS over HypoARDS,
# average log2 CPM, moderated t, p value and Benjamini-Hochberg FDR.

library(tidyverse)
library(edgeR)

clean_directory <- file.path(here::here(), "data", "geo", "clean")
results_directory <- file.path(here::here(), "data", "results", "ajrccm_2025")
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

model_subjects <-
  read_csv(file.path(clean_directory, "noncovid_bulk_subjects.csv"),
           show_col_types = FALSE) %>%
  mutate(phenotype = factor(phenotype, levels = c("HypoARDS", "HyperARDS")))
gene_annotation <- read_csv(file.path(clean_directory, "gene_annotation.csv"),
                            show_col_types = FALSE)
counts <-
  read_csv(file.path(clean_directory, "noncovid_bulk_counts.csv"), show_col_types = FALSE) %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix()
stopifnot(identical(colnames(counts), model_subjects$sample_id))

if (anyNA(model_subjects$age)) {
  stop("Age is missing for ", sum(is.na(model_subjects$age)), " samples. The models adjust ",
       "for age, which must be available from GEO; see script 0.")
}

# HypoARDS is the reference level, so the phenotype coefficient is the log2
# fold change of HyperARDS over HypoARDS.
design <- model.matrix(~ rna_batch + phenotype + sex + age, data = model_subjects)

# Normalization factors are computed on all non-zero genes before filtering,
# and filtering keeps the original library sizes, as in rnaseq_to_limma().
expression_list <- normLibSizes(DGEList(counts, remove.zeros = TRUE))
expressed_genes <- filterByExpr(expression_list, design = model.matrix(~ 1, model_subjects))

phenotype_fit <- voomLmFit(expression_list[expressed_genes, ], design = design) %>%
  eBayes(robust = TRUE)

phenotype_results <-
  topTable(phenotype_fit, coef = "phenotypeHyperARDS", number = Inf, sort.by = "P") %>%
  rownames_to_column("ensembl_gene_id") %>%
  as_tibble() %>%
  left_join(gene_annotation, by = "ensembl_gene_id") %>%
  transmute(ensembl_gene_id, gene_symbol, log2_fold_change = logFC,
            average_log2_cpm = AveExpr, moderated_t = t, p_value = P.Value,
            fdr = adj.P.Val)

write_csv(phenotype_results, file.path(results_directory, "noncovid_bulk_hyper_vs_hypo.csv"))

message("Non-COVID bulk RNA-seq: ", nrow(model_subjects), " samples (",
        paste(names(table(model_subjects$phenotype)), table(model_subjects$phenotype),
              collapse = ", "), "); ", sum(expressed_genes), " of ", nrow(counts),
        " genes passed the expression filter")
message("Genes with FDR < 0.05: ", sum(phenotype_results$fdr < 0.05),
        "; FDR < 0.1: ", sum(phenotype_results$fdr < 0.1))
