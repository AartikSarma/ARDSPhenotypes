# ARDS Phenotypes - 4 - Reproducible genes and gene set enrichment.R
#
# Inputs are the differential expression results reported in the paper:
#   bulk: published_results/Supplementary Data E1 - Bulk RNASeq.xlsx. The
#     published models used exact ages and voom weights from samples not on
#     GEO, so this table is used rather than the refits of scripts 1 and 2;
#   single-cell: Supplementary Data E5, written by script 3 to
#     data/results/ajrccm_2025/.
# Fold changes are HyperARDS relative to HypoARDS throughout.
#
# 1. Reproducible genes: Ensembl genes with FDR < 0.1 in both bulk cohorts. Each is
#    called up in HyperARDS or in HypoARDS by the sign of its non-COVID fold
#    change, and described from MyGene.info (Supplementary Data E2).
# 2. Reactome gene set enrichment (fgseaMultilevel) on genes ranked by log2
#    fold change, with sets of at least 30 genes (tidyseq v0.1.2
#    dge_to_gsea()), run the same way for bulk and single-cell results
#    (Supplementary Data E4). Bulk fold changes are averaged per gene symbol;
#    single-cell sheets report only pathways with FDR < 0.1. The published
#    single-cell analysis used sets of at least 15 genes.
#    Reactome sets come from msigdbr 7.5.1 (MSigDB 7.5.1), the release used in
#    2024; newer releases will change the sets and the results.
#
#    fgsea estimates p values by random sampling, so every call is seeded,
#    with many permutations (nPermSimple = 100000) and a precise multilevel
#    estimate (sampleSize = 1001, eps = 0), because several pathways lie
#    close to FDR 0.1.

# 3. Pathways for Figures 2B and E4: pathways with FDR < 0.1 and the same sign
#    in both bulk cohorts, joined to the single-cell scores.
#
# Run from the repository root after script 3. Writes to
# data/results/ajrccm_2025/:
#   Supplementary Data E2 - Gene descriptions.xlsx
#   Supplementary Data E4 - GSEA results.xlsx
#   reproducible_genes.csv                    the genes behind Figure 1C
#   bulk_gsea_shared_pathways.csv             Figure 2B data
#   bulk_and_single_cell_gsea.csv             Figure E4 data

library(tidyverse)
library(readxl)
library(fgsea)
library(msigdbr)

published_directory <- file.path(here::here(), "published_results")
results_directory <- file.path(here::here(), "data", "results", "ajrccm_2025")
dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)

bulk_path <- file.path(published_directory, "Supplementary Data E1 - Bulk RNASeq.xlsx")
single_cell_path <- file.path(results_directory, "Supplementary Data E5 - scSeq MAST results.xlsx")

fdr_threshold <- 0.1

stopifnot(packageVersion("msigdbr") == "7.5.1")

# Published differential expression ----------------------------------------

bulk_results <- list(
  "non-COVID" = read_excel(bulk_path, sheet = "S1A - ALI Limma"),
  "COVID" = read_excel(bulk_path, sheet = "S1B - COMET Dream")
) %>%
  map(~ transmute(.x, ensembl_gene_id = ensembl, gene_symbol = symbol,
                  log2_fold_change = logFC, p_value = `P-value`, fdr = FDR))

# Supplementary Data E5 sheet for each cell type, in the order of Data E4.
single_cell_sheets <- c(
  "Alveolar mac" = "S5C - Alveolar macrophages",
  "Monocyte" = "S5D - Monocytes",
  "Monocyte-derived mac" = "S5B - Monocyte-derived macs",
  "Neutrophil" = "S5E - Neutrophils",
  "T cell" = "S5A - T cells"
)
single_cell_results <- map(single_cell_sheets, function(sheet) {
  read_excel(single_cell_path, sheet = sheet) %>%
    transmute(gene_symbol = HGNC, log2_fold_change = `Log2 Fold Change`,
              p_value = `P-value`, fdr = FDR)
})

# Reproducible genes --------------------------------------------------------

# Genes are matched by Ensembl ID, as some gene symbols name more than one
# Ensembl gene.
significant_genes <- map(bulk_results, ~ .x$ensembl_gene_id[.x$fdr < fdr_threshold])
reproducible_gene_ids <- reduce(significant_genes, intersect)

reproducible_genes <-
  bulk_results[["non-COVID"]] %>%
  filter(ensembl_gene_id %in% reproducible_gene_ids) %>%
  dplyr::select(ensembl_gene_id, gene_symbol, noncovid_log2_fold_change = log2_fold_change,
                noncovid_fdr = fdr) %>%
  left_join(bulk_results[["COVID"]] %>%
              dplyr::select(ensembl_gene_id, covid_log2_fold_change = log2_fold_change,
                            covid_fdr = fdr),
            by = "ensembl_gene_id") %>%
  mutate(direction = if_else(noncovid_log2_fold_change > 0, "Up in HyperARDS", "Up in HypoARDS"))

write_csv(reproducible_genes, file.path(results_directory, "reproducible_genes.csv"))

# Gene summaries from MyGene.info, which change as the annotation is updated.
gene_summaries <-
  mygene::getGenes(reproducible_genes$ensembl_gene_id, fields = "summary") %>%
  as.data.frame() %>%
  as_tibble() %>%
  transmute(ensembl_gene_id = query, description = summary)

gene_descriptions <-
  reproducible_genes %>%
  left_join(gene_summaries, by = "ensembl_gene_id") %>%
  transmute(direction, ensembl = ensembl_gene_id, symbol = gene_symbol, description)

writexl::write_xlsx(
  list("S2A - Up in HyperARDS" = gene_descriptions %>%
         filter(direction == "Up in HyperARDS") %>% dplyr::select(-direction),
       "S2B - Up in HypoARDS" = gene_descriptions %>%
         filter(direction == "Up in HypoARDS") %>% dplyr::select(-direction)),
  file.path(results_directory, "Supplementary Data E2 - Gene descriptions.xlsx"))

# Gene set enrichment ------------------------------------------------------

reactome_gene_sets <-
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  { split(.$gene_symbol, .$gs_name) }

clean_pathway_name <- function(pathway) {
  pathway %>% str_remove("^REACTOME_") %>% str_replace_all("_", " ")
}

run_gsea <- function(results) {
  ranked_fold_changes <-
    results %>%
    filter(!is.na(log2_fold_change), !is.na(gene_symbol)) %>%
    group_by(gene_symbol) %>%
    summarize(log2_fold_change = mean(log2_fold_change)) %>%
    deframe()
  set.seed(20240724)
  fgseaMultilevel(pathways = reactome_gene_sets, stats = ranked_fold_changes,
                  minSize = 30, nPermSimple = 100000, sampleSize = 1001, eps = 0) %>%
    as_tibble()
}

bulk_gsea <- map(bulk_results, run_gsea)
single_cell_gsea <- map(single_cell_results, ~ run_gsea(.x) %>% filter(padj < fdr_threshold))

tidy_gsea <- function(gsea_results) {
  gsea_results %>%
    transmute(pathway = clean_pathway_name(pathway), NES, `P-value` = pval, fdr = padj,
              leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ", ")))
}

gsea_sheets <- list(
  "S4A - ALI BulkRNA" = bulk_gsea[["non-COVID"]],
  "S4B - COMET BulkRNA" = bulk_gsea[["COVID"]],
  "S4C - Alv Macrophage" = single_cell_gsea[["Alveolar mac"]],
  "S4D - Monocytes" = single_cell_gsea[["Monocyte"]],
  "S4E - Monocyte-derived mac" = single_cell_gsea[["Monocyte-derived mac"]],
  "S4F - Neutrophils" = single_cell_gsea[["Neutrophil"]],
  "S4G - T cells" = single_cell_gsea[["T cell"]]
) %>%
  map(tidy_gsea)

writexl::write_xlsx(gsea_sheets,
                    file.path(results_directory, "Supplementary Data E4 - GSEA results.xlsx"))

# Pathways for Figures 2B and E4 ------------------------------------------

bulk_gsea_shared_pathways <-
  gsea_sheets[["S4A - ALI BulkRNA"]] %>%
  dplyr::select(pathway, noncovid_nes = NES, noncovid_fdr = fdr) %>%
  inner_join(gsea_sheets[["S4B - COMET BulkRNA"]] %>%
               dplyr::select(pathway, covid_nes = NES, covid_fdr = fdr),
             by = "pathway") %>%
  filter(noncovid_fdr < fdr_threshold, covid_fdr < fdr_threshold,
         sign(noncovid_nes) == sign(covid_nes)) %>%
  arrange(desc(noncovid_nes * covid_nes * sign(noncovid_nes)))

write_csv(bulk_gsea_shared_pathways,
          file.path(results_directory, "bulk_gsea_shared_pathways.csv"))

# Pathways shared by the bulk cohorts that are also enriched, in the same
# direction as in the non-COVID cohort, in at least one cell type.
single_cell_nes <-
  imap(single_cell_gsea, ~ tidy_gsea(.x) %>% transmute(pathway, comparison = .y, NES)) %>%
  list_rbind()

bulk_and_single_cell_gsea <-
  bulk_gsea_shared_pathways %>%
  transmute(pathway, `non-COVID` = noncovid_nes, COVID = covid_nes) %>%
  pivot_longer(-pathway, names_to = "comparison", values_to = "NES") %>%
  bind_rows(semi_join(single_cell_nes, bulk_gsea_shared_pathways, by = "pathway")) %>%
  group_by(pathway) %>%
  filter(any(comparison %in% names(single_cell_gsea) &
               sign(NES) == sign(NES[comparison == "non-COVID"]))) %>%
  ungroup() %>%
  mutate(sequencing = if_else(comparison %in% c("non-COVID", "COVID"),
                              "Bulk RNA-seq", "Single-cell RNA-seq"))

write_csv(bulk_and_single_cell_gsea,
          file.path(results_directory, "bulk_and_single_cell_gsea.csv"))

message(nrow(reproducible_genes), " reproducible genes (",
        sum(reproducible_genes$direction == "Up in HyperARDS"), " up in HyperARDS, ",
        sum(reproducible_genes$direction == "Up in HypoARDS"), " up in HypoARDS)")
message(nrow(bulk_gsea_shared_pathways), " Reactome pathways enriched in the same direction ",
        "in both bulk cohorts (FDR < ", fdr_threshold, ")")
