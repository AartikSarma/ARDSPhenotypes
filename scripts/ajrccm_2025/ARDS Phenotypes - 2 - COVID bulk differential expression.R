# ARDS Phenotypes - 2 - COVID bulk differential expression.R
#
# Tests for genes differentially expressed between HyperARDS and HypoARDS in
# tracheal aspirate bulk RNA-seq from the COVID cohort (GSE350140). Subjects
# contribute samples from more than one ventilator day, so the model is a
# linear mixed model fit with dream (variancePartition), with dream's default
# ("adaptive") degrees of freedom, which are Satterthwaite's for 12 or more
# samples:
#
#   ~ phenotype + scale(age) + sex + (1 | rna_batch/subject_id)
#
# Subjects are nested within RNA-seq batch. Three subjects have samples in two
# batches; the nesting gives each of them a separate intercept per batch, as in
# the published fit.
#
# Genes with zero counts in every sample are removed, library sizes are
# upper-quartile normalized, genes are kept by filterByExpr() with no
# grouping, and dream's fit is moderated with (non-robust) empirical Bayes.
# The published analysis filtered genes and estimated voom weights before
# removing PCA outlier samples; GEO holds only the 72 samples that remained,
# so here both steps use those 72.
#
# Age is top-coded at 90 (script 0) and scaled across samples, not subjects,
# as in the published fit.
#
# Because voom weights and the gene filter are computed without the outlier
# samples, and because age is top-coded, results may differ subtly from those
# published. Additional data are available on reasonable request.
#
# Run from the repository root after script 0. Writes
#   data/results/ajrccm_2025/covid_bulk_hyper_vs_hypo.csv
# with one row per tested gene: log2 fold change of HyperARDS over HypoARDS,
# average log2 CPM, moderated t, its degrees of freedom, p value and
# Benjamini-Hochberg FDR.

library(tidyverse)
library(edgeR)
library(variancePartition)
library(BiocParallel)

clean_directory <- file.path(here::here(), "data", "geo", "clean")
results_directory <- file.path(here::here(), "data", "results", "ajrccm_2025")
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

model_subjects <-
  read_csv(file.path(clean_directory, "covid_bulk_subjects.csv"),
           show_col_types = FALSE) %>%
  mutate(phenotype = factor(phenotype, levels = c("HypoARDS", "HyperARDS")),
         scale_age = as.numeric(scale(age))) %>%
  as.data.frame()
gene_annotation <- read_csv(file.path(clean_directory, "gene_annotation.csv"),
                            show_col_types = FALSE)
counts <-
  read_csv(file.path(clean_directory, "covid_bulk_counts.csv"), show_col_types = FALSE) %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix()
stopifnot(identical(colnames(counts), model_subjects$sample_id))
rownames(model_subjects) <- model_subjects$sample_id

parallel_workers <- MulticoreParam(max(1, parallel::detectCores() - 2))

# HypoARDS is the reference level, so the phenotype coefficient is the log2
# fold change of HyperARDS over HypoARDS.
model_formula <- ~ phenotype + scale_age + sex + (1 | rna_batch/subject_id)

expression_list <- normLibSizes(DGEList(counts, remove.zeros = TRUE), method = "upperquartile")
expressed_genes <- filterByExpr(expression_list)
voom_expression <- voomWithDreamWeights(expression_list[expressed_genes, ], model_formula,
                                        model_subjects,
                                        BPPARAM = parallel_workers)

phenotype_fit <- dream(voom_expression, model_formula, model_subjects,
                       ddf = "adaptive", BPPARAM = parallel_workers)
failed_model_errors <- attr(phenotype_fit, "errors")
phenotype_fit <- eBayes(phenotype_fit)
phenotype_degrees_of_freedom <- phenotype_fit$df.total[, "phenotypeHyperARDS"]

phenotype_results <-
  topTable(phenotype_fit, coef = "phenotypeHyperARDS", number = Inf, sort.by = "P") %>%
  rownames_to_column("ensembl_gene_id") %>%
  as_tibble() %>%
  left_join(gene_annotation, by = "ensembl_gene_id") %>%
  transmute(ensembl_gene_id, gene_symbol, log2_fold_change = logFC,
            average_log2_cpm = AveExpr, moderated_t = t,
            degrees_of_freedom = phenotype_degrees_of_freedom[ensembl_gene_id],
            p_value = P.Value, fdr = adj.P.Val)

write_csv(phenotype_results, file.path(results_directory, "covid_bulk_hyper_vs_hypo.csv"))

message("COVID bulk RNA-seq: ", nrow(model_subjects), " samples from ",
        n_distinct(model_subjects$subject_id), " subjects; ", sum(expressed_genes), " of ",
        nrow(counts), " genes passed the expression filter; the mixed model failed for ",
        length(failed_model_errors), " of them")
if (length(failed_model_errors) > 0) {
  message("Model errors, by type:")
  unlist(failed_model_errors) %>%
    str_remove_all("-?[0-9.]+e[-+][0-9]+|\\s+$") %>%
    table() %>%
    print()
}
message("Genes with FDR < 0.1: ", sum(phenotype_results$fdr < 0.1))
