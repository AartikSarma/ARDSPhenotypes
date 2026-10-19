# ARDS Phenotypes - 0 - Import and clean GEO data.R
#
# Downloads the public data deposited in GEO (SuperSeries GSE200849) and writes
# cleaned, analysis-ready tables. The later scripts in this folder, which
# reproduce the AJRCCM 2025 paper, start from these outputs.
#
#   GSE200847  non-COVID ARDS, tracheal aspirate bulk RNA-seq (41 samples)
#   GSE350140  COVID ARDS, tracheal aspirate bulk RNA-seq (72 samples from
#              repeated ventilator days)
#   GSE200848  ARDS, tracheal aspirate single-cell RNA-seq (one Seurat object)
#
# Run from the repository root. Inputs are downloaded to data/geo/raw/;
# outputs are written to data/geo/clean/:
#   gene_annotation.csv            Ensembl v99 gene symbol and biotype
#   noncovid_bulk_subjects.csv     one row per sample, GSE200847
#   noncovid_bulk_counts.csv       genes x samples, columns ordered as subjects
#   covid_bulk_subjects.csv        one row per sample, GSE350140
#   covid_bulk_counts.csv          genes x samples, columns ordered as subjects
#   scseq_subjects.csv             one row per GEO sample record, GSE200848
#   scseq_seurat.rds               Seurat object, updated to the installed Seurat
#
# Age, a covariate of the bulk models, is taken from GEO, which records it as
# a whole number of years or as ">89" for patients older than 89 (HIPAA Safe
# Harbor). ">89" and any age above 89 become 90, flagged by age_top_coded.

library(tidyverse)
library(GEOquery)
library(Seurat)

project_directory <- here::here()
raw_directory <- file.path(project_directory, "data", "geo", "raw")
clean_directory <- file.path(project_directory, "data", "geo", "clean")
dir.create(raw_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(clean_directory, recursive = TRUE, showWarnings = FALSE)

options(timeout = 3600) # the single-cell object is 1.6 GB

# Download -----------------------------------------------------------------

download_if_missing <- function(url, destination_path) {
  if (!file.exists(destination_path)) {
    message("Downloading ", basename(destination_path))
    download.file(url, destination_path, mode = "wb")
  }
  destination_path
}

geo_supplementary_url <- function(accession, file_name) {
  series_stub <- str_replace(accession, "\\d{3}$", "nnn")
  paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", series_stub, "/",
         accession, "/suppl/", file_name)
}

# getGEO reuses a series matrix already present in destdir.
read_geo_sample_annotations <- function(accession) {
  expression_set <- getGEO(accession, destdir = raw_directory,
                           GSEMatrix = TRUE, getGPL = FALSE)[[1]]
  pData(expression_set) %>%
    as_tibble() %>%
    rename_with(~ .x %>% str_remove(":ch1$") %>% str_to_lower() %>%
                  str_replace_all("[^a-z0-9]+", "_"))
}

noncovid_counts_path <- download_if_missing(
  geo_supplementary_url("GSE200847", "GSE200847_cohort1.counts.csv.gz"),
  file.path(raw_directory, "GSE200847_cohort1.counts.csv.gz"))
covid_counts_path <- download_if_missing(
  geo_supplementary_url("GSE350140", "GSE350140_cohort2.counts.csv.gz"),
  file.path(raw_directory, "GSE350140_cohort2.counts.csv.gz"))
scseq_path <- download_if_missing(
  geo_supplementary_url("GSE200848", "GSE200848_cohort3.counts.Rds"),
  file.path(raw_directory, "GSE200848_cohort3.counts.Rds"))

# The data processing notes name Ensembl v99 as the pseudo-alignment
# reference, so symbols come from that release rather than current BioMart.
gtf_path <- download_if_missing(
  "https://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz",
  file.path(raw_directory, "Homo_sapiens.GRCh38.99.gtf.gz"))

noncovid_annotations <- read_geo_sample_annotations("GSE200847")
covid_annotations <- read_geo_sample_annotations("GSE350140")
scseq_annotations <- read_geo_sample_annotations("GSE200848")

# Gene annotation ----------------------------------------------------------

# gene_symbol is the Ensembl gene_name: the HGNC symbol where one exists,
# otherwise a clone-based name.
gene_annotation <-
  read_tsv(gtf_path, comment = "#", col_names = FALSE,
           col_types = cols(.default = "c")) %>%
  filter(X3 == "gene") %>%
  transmute(
    ensembl_gene_id = str_match(X9, 'gene_id "([^"]+)"')[, 2],
    gene_symbol = str_match(X9, 'gene_name "([^"]+)"')[, 2],
    gene_biotype = str_match(X9, 'gene_biotype "([^"]+)"')[, 2]
  )

# Bulk RNA-seq -------------------------------------------------------------

# Count tables carry samples in a different order from the series matrix, so
# columns are matched by name. Counts are tximport scaledTPM values.
read_bulk_counts <- function(counts_path, sample_ids) {
  counts <- read_csv(counts_path, show_col_types = FALSE)
  stopifnot(
    setequal(setdiff(colnames(counts), "ensembl_gene_id"), sample_ids),
    !anyDuplicated(counts$ensembl_gene_id),
    all(counts$ensembl_gene_id %in% gene_annotation$ensembl_gene_id)
  )
  count_values <- as.matrix(counts[, -1])
  stopifnot(!anyNA(count_values), all(count_values >= 0),
            all(count_values == round(count_values)))
  counts %>% dplyr::select(ensembl_gene_id, all_of(sample_ids))
}

phenotype_levels <- c("HypoARDS", "HyperARDS")

top_coded_age <- 90
parse_geo_age <- function(geo_age) {
  geo_age <- str_trim(as.character(geo_age))
  stopifnot(all(geo_age == ">89" | str_detect(geo_age, "^\\d+$")))
  numeric_age <- if_else(geo_age == ">89", top_coded_age, suppressWarnings(as.numeric(geo_age)))
  pmin(numeric_age, top_coded_age)
}

noncovid_bulk_subjects <-
  noncovid_annotations %>%
  transmute(
    sample_id = title,
    subject_id = title,
    geo_accession,
    phenotype = factor(phenotype, levels = phenotype_levels),
    sex,
    age = parse_geo_age(age),
    age_top_coded = age == top_coded_age,
    rna_batch = batch,
    sofa_score = parse_integer(sofa, na = "NA"),
    vasopressor = case_when(vasopressor == "Yes" ~ TRUE, vasopressor == "No" ~ FALSE)
  ) %>%
  arrange(sample_id)

noncovid_bulk_counts <- read_bulk_counts(noncovid_counts_path,
                                         noncovid_bulk_subjects$sample_id)

# Titles encode subject and ventilator day: COMET1145_VENTDAY2.
# Several subjects contribute samples from more than one day.
covid_bulk_subjects <-
  covid_annotations %>%
  transmute(
    sample_id = title,
    subject_id = str_extract(title, "^[^_]+"),
    vent_day = parse_integer(str_extract(title, "(?<=VENTDAY)\\d+$")),
    geo_accession,
    phenotype = factor(phenotype, levels = phenotype_levels),
    sex,
    age = parse_geo_age(age),
    age_top_coded = age == top_coded_age,
    rna_batch = batch,
    sofa_score = parse_integer(sofa, na = "NA"),
    vasopressor = case_when(vasopressor == "1" ~ TRUE, vasopressor == "0" ~ FALSE)
  ) %>%
  arrange(subject_id, vent_day)

# Phenotype and age are subject characteristics, so they must not vary across
# days.
stopifnot(covid_bulk_subjects %>% distinct(subject_id, age) %>%
            pull(subject_id) %>% anyDuplicated() == 0)
stopifnot(!anyNA(covid_bulk_subjects$vent_day),
          !anyDuplicated(covid_bulk_subjects$sample_id),
          covid_bulk_subjects %>% distinct(subject_id, phenotype) %>%
            pull(subject_id) %>% anyDuplicated() == 0)

covid_bulk_counts <- read_bulk_counts(covid_counts_path,
                                      covid_bulk_subjects$sample_id)

# Missing vasopressor and SOFA values must come from a literal "NA" in GEO,
# not from an unrecognized value.
stopifnot(
  !anyNA(noncovid_bulk_subjects$phenotype), !anyNA(covid_bulk_subjects$phenotype),
  sum(is.na(noncovid_bulk_subjects$vasopressor)) == sum(noncovid_annotations$vasopressor == "NA"),
  sum(is.na(covid_bulk_subjects$vasopressor)) == sum(covid_annotations$vasopressor == "NA"),
  sum(is.na(covid_bulk_subjects$sofa_score)) == sum(covid_annotations$sofa == "NA")
)

# Single-cell RNA-seq ------------------------------------------------------

# The single-cell subjects are a separate cohort from the non-COVID bulk
# subjects.
scseq_subjects <-
  scseq_annotations %>%
  transmute(
    sample_id = title,
    geo_accession,
    phenotype = factor(phenotype, levels = phenotype_levels),
    sex,
    vent_day = parse_integer(vent_day),
    sofa_score = parse_integer(sofa, na = "NA"),
    vasopressor = case_when(vasopressor == "Yes" ~ TRUE, vasopressor == "No" ~ FALSE)
  )

scseq_seurat <- UpdateSeuratObject(readRDS(scseq_path))

# Library names appear as both HS* and sc*; key cells to the GEO titles (HS*).
# The object must hold exactly the subjects in the GEO sample records, and its
# own phenotype labels must agree with them; phenotype is then taken from the
# records.
scseq_seurat$sample_id <- str_replace(scseq_seurat$orig.ident, "^sc", "HS")
stopifnot(setequal(unique(scseq_seurat$sample_id), scseq_subjects$sample_id),
          scseq_seurat[[]] %>%
            left_join(scseq_subjects, by = "sample_id",
                      suffix = c("_object", "_geo")) %>%
            with(all(phenotype_object == as.character(phenotype_geo))))
scseq_seurat$phenotype <- as.character(
  scseq_subjects$phenotype[match(scseq_seurat$sample_id, scseq_subjects$sample_id)])

scseq_subjects <-
  scseq_subjects %>%
  left_join(scseq_seurat[[]] %>% count(sample_id, name = "cells"), by = "sample_id")

# Write --------------------------------------------------------------------

write_csv(gene_annotation, file.path(clean_directory, "gene_annotation.csv"))
write_csv(noncovid_bulk_subjects, file.path(clean_directory, "noncovid_bulk_subjects.csv"))
write_csv(noncovid_bulk_counts, file.path(clean_directory, "noncovid_bulk_counts.csv"))
write_csv(covid_bulk_subjects, file.path(clean_directory, "covid_bulk_subjects.csv"))
write_csv(covid_bulk_counts, file.path(clean_directory, "covid_bulk_counts.csv"))
write_csv(scseq_subjects, file.path(clean_directory, "scseq_subjects.csv"))
saveRDS(scseq_seurat, file.path(clean_directory, "scseq_seurat.rds"))

# Summary ------------------------------------------------------------------

summarize_bulk_cohort <- function(bulk_subjects, cohort_name) {
  message(cohort_name, ": ", nrow(bulk_subjects), " samples from ",
          n_distinct(bulk_subjects$subject_id), " subjects; age top-coded at ", top_coded_age,
          " for ", n_distinct(bulk_subjects$subject_id[bulk_subjects$age_top_coded]), " subjects")
  bulk_subjects %>%
    distinct(subject_id, phenotype) %>%
    count(phenotype, name = "subjects") %>%
    print()
}

summarize_bulk_cohort(noncovid_bulk_subjects, "Non-COVID bulk RNA-seq (GSE200847)")
summarize_bulk_cohort(covid_bulk_subjects, "COVID bulk RNA-seq (GSE350140)")

message("Biotypes of the genes in the bulk count tables:")
gene_annotation %>%
  semi_join(noncovid_bulk_counts, by = "ensembl_gene_id") %>%
  count(gene_biotype, name = "genes") %>%
  print()

message("Single-cell RNA-seq (GSE200848): ", ncol(scseq_seurat), " cells from ",
        n_distinct(scseq_seurat$sample_id), " subjects")
scseq_subjects %>% dplyr::select(sample_id, phenotype, cells) %>% print()
# The deposited object holds only cells that passed quality control
# (300 < nCount_RNA < 30000, percent.mt < 10).
message("nCount_RNA range: ", paste(range(scseq_seurat$nCount_RNA), collapse = " to "),
        "; maximum percent.mt: ", round(max(scseq_seurat$percent.mt), 2))
