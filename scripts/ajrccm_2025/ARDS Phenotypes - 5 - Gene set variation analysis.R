# ARDS Phenotypes - 5 - Gene set variation analysis.R
#
# Scores every bulk RNA-seq sample for a HyperARDS signature derived from the
# non-COVID cohort and tests whether the score differs between phenotypes in
# both cohorts (Figures 2A and E2, Supplementary Data E3).
#
# The signature is the 100 genes with the smallest p values among those up in
# HyperARDS in the non-COVID cohort (log2 fold change > 0.5, FDR < 0.1), taken
# from the published Supplementary Data E1. Samples are scored by GSVA on voom
# log2 CPM values, normalized by upper-quartile scaling:
#   - non-COVID: genes kept by filterByExpr() before normalization;
#   - COVID: genes restricted to the 13,024 published genes (E1, sheet S1B),
#     so that GSVA ranks the same genes as the published voom object.
# The published COVID scores came from that voom object, which also held the
# four PCA outlier samples that GEO lacks. GSVA estimates each gene's
# distribution across samples, so COVID scores differ subtly from those
# published; additional data are available on reasonable request.
#
# Phenotype is tested with the covariates of the differential expression
# models, age entering unscaled:
#   non-COVID:  lm(score ~ phenotype + age + sex + rna_batch)
#   COVID:      lmer(score ~ phenotype + age + sex + (1 | rna_batch/subject_id)),
#               Satterthwaite degrees of freedom (lmerTest)
#   COVID, first sample per subject (Figure E2):
#               lmer(score ~ phenotype + age + sex + (1 | rna_batch))
# Age comes from GEO (script 0) and is top-coded at 90.
#
# Run from the repository root after script 0. Writes to
# data/results/ajrccm_2025/:
#   Supplementary Data E3 - GSVA.xlsx
#   gsva_phenotype_tests.csv

library(tidyverse)
library(readxl)
library(edgeR)
library(GSVA)
library(lmerTest)

clean_directory <- file.path(here::here(), "data", "geo", "clean")
published_directory <- file.path(here::here(), "published_results")
results_directory <- file.path(here::here(), "data", "results", "ajrccm_2025")
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

bulk_path <- file.path(published_directory, "Supplementary Data E1 - Bulk RNASeq.xlsx")

read_counts <- function(file_name) {
  read_csv(file.path(clean_directory, file_name), show_col_types = FALSE) %>%
    column_to_rownames("ensembl_gene_id") %>%
    as.matrix()
}

noncovid_subjects <-
  read_csv(file.path(clean_directory, "noncovid_bulk_subjects.csv"), show_col_types = FALSE) %>%
  mutate(phenotype = factor(phenotype, levels = c("HypoARDS", "HyperARDS")))
covid_subjects <-
  read_csv(file.path(clean_directory, "covid_bulk_subjects.csv"), show_col_types = FALSE) %>%
  mutate(phenotype = factor(phenotype, levels = c("HypoARDS", "HyperARDS")))
noncovid_counts <- read_counts("noncovid_bulk_counts.csv")
covid_counts <- read_counts("covid_bulk_counts.csv")
stopifnot(identical(colnames(noncovid_counts), noncovid_subjects$sample_id),
          identical(colnames(covid_counts), covid_subjects$sample_id))

# Signature -----------------------------------------------------------------

noncovid_published <- read_excel(bulk_path, sheet = "S1A - ALI Limma")
covid_published_genes <- read_excel(bulk_path, sheet = "S1B - COMET Dream")$ensembl

hyperards_signature <- list(
  hyperards_signature =
    noncovid_published %>%
    filter(logFC > 0.5, FDR < 0.1) %>%
    arrange(`P-value`) %>%
    head(100) %>%
    pull(ensembl)
)

# Scores --------------------------------------------------------------------

noncovid_expression <- DGEList(noncovid_counts, remove.zeros = TRUE)
noncovid_expression <-
  noncovid_expression[filterByExpr(noncovid_expression), ] %>%
  normLibSizes(method = "upperquartile") %>%
  voom()

covid_expression <-
  DGEList(covid_counts, remove.zeros = TRUE) %>%
  normLibSizes(method = "upperquartile")
covid_expression <- voom(covid_expression[covid_published_genes, ])

score_samples <- function(voom_expression) {
  gsva(gsvaParam(voom_expression$E, hyperards_signature), verbose = FALSE)[1, ] %>%
    enframe("sample_id", "gsva_score")
}

noncovid_scores <- noncovid_subjects %>%
  left_join(score_samples(noncovid_expression), by = "sample_id")
covid_scores <- covid_subjects %>%
  left_join(score_samples(covid_expression), by = "sample_id")

writexl::write_xlsx(
  list("S3A - ALI" = noncovid_scores %>%
         transmute(subjectid = subject_id, phenotype = as.character(phenotype),
                   `GSVA score` = gsva_score),
       "S3B - COMET" = covid_scores %>%
         transmute(subject = subject_id, phenotype = as.character(phenotype), vent_day,
                   `GSVA score` = gsva_score)),
  file.path(results_directory, "Supplementary Data E3 - GSVA.xlsx"))

# Phenotype tests ----------------------------------------------------------

phenotype_term <- function(model, cohort) {
  summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    filter(term == "phenotypeHyperARDS") %>%
    transmute(cohort, difference = Estimate, standard_error = `Std. Error`,
              p_value = .[[ncol(.)]])
}

covid_first_samples <- covid_scores %>%
  arrange(subject_id, vent_day) %>%
  distinct(subject_id, .keep_all = TRUE)

gsva_phenotype_tests <- bind_rows(
  lm(gsva_score ~ phenotype + age + sex + rna_batch, data = noncovid_scores) %>%
    phenotype_term("non-COVID"),
  lmer(gsva_score ~ phenotype + age + sex + (1 | rna_batch/subject_id), data = covid_scores) %>%
    phenotype_term("COVID, all samples"),
  lmer(gsva_score ~ phenotype + age + sex + (1 | rna_batch), data = covid_first_samples) %>%
    phenotype_term("COVID, first sample per subject")
)

write_csv(gsva_phenotype_tests, file.path(results_directory, "gsva_phenotype_tests.csv"))

message("HyperARDS signature: ", length(hyperards_signature$hyperards_signature), " genes")
message("Phenotype tests of the GSVA score (HyperARDS minus HypoARDS):")
print(gsva_phenotype_tests, digits = 3)
