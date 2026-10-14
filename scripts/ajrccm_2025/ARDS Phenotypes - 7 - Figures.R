# ARDS Phenotypes - 7 - Figures.R
#
# Draws the transcriptomic figures of the AJRCCM 2025 paper:
#   Figure 1A, 1B  volcano plots of bulk differential expression (Data E1)
#   Figure 1C      heatmaps of the reproducible genes in each bulk cohort
#   (unnumbered)   expression of the reproducible genes by cell type and
#                  phenotype in the single-cell RNA-seq
#   (unnumbered)   UMAPs of the single-cell RNA-seq showing the expression of
#                  each reproducible gene up in HyperARDS (PNG, as the plot
#                  holds about 26,000 points per panel)
#   (unnumbered)   per-cell scores for those genes (all, ribosomal, other), by
#                  phenotype: on the UMAP, over all cells and within each cell
#                  type, and against each subject's T and NK cell fraction
#   Figure 2A      HyperARDS signature GSVA scores by phenotype
#   Figure 2B      Reactome pathways enriched in the same direction in both
#                  bulk cohorts (a curated subset)
#   Figure E2      GSVA scores, first COVID sample per subject
#   Figure E3      the shared bulk pathways that are also enriched in at least
#                  one single-cell type
#
# Phenotypes are coloured from the Okabe-Ito palette. Continuous values that
# run from low to high use viridis palettes; values centred on zero (Figure E3
# NES, scaled expression and the HyperARDS gene score) use diverging scales.
# Heatmap dendrograms are ordered by dendsort.
#
# Run from the repository root after scripts 0, 3, 4 and 5. Writes PDFs to
# data/results/ajrccm_2025/figures/.

library(tidyverse)
library(readxl)
library(DESeq2)
library(pheatmap)
library(dendsort)
library(ggrepel)
library(patchwork)
library(viridis)
library(SeuratObject)

clean_directory <- file.path(here::here(), "data", "geo", "clean")
published_directory <- file.path(here::here(), "published_results")
results_directory <- file.path(here::here(), "data", "results", "ajrccm_2025")
figure_directory <- file.path(results_directory, "figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

# Okabe-Ito vermillion and blue.
phenotype_colours <- c(HyperARDS = "#D55E00", HypoARDS = "#0072B2")
phenotype_labels <- c(HypoARDS = "Hypoinflammatory", HyperARDS = "Hyperinflammatory")
cohort_titles <- c("non-COVID" = "ALI - Non-COVID ARDS", "COVID" = "COMET - COVID ARDS")

fdr_threshold <- 0.1
fold_change_threshold <- 0.5

bulk_path <- file.path(published_directory, "Supplementary Data E1 - Bulk RNASeq.xlsx")
bulk_results <- list(
  "non-COVID" = read_excel(bulk_path, sheet = "S1A - ALI Limma"),
  "COVID" = read_excel(bulk_path, sheet = "S1B - COMET Dream")
)

# Figure 1A, 1B: volcano plots ---------------------------------------------

# Genes with FDR < 0.1 and an absolute log2 fold change above 0.5 are coloured
# by direction; the ten with the smallest p values in each direction are
# labelled.
plot_volcano <- function(results, title) {
  results <- results %>%
    mutate(direction = case_when(
      FDR < fdr_threshold & logFC > fold_change_threshold ~ "HyperARDS",
      FDR < fdr_threshold & logFC < -fold_change_threshold ~ "HypoARDS",
      TRUE ~ "Not significant"))
  labelled_genes <- results %>%
    filter(direction != "Not significant", !is.na(symbol)) %>%
    arrange(`P-value`) %>%
    group_by(direction) %>%
    slice_head(n = 10) %>%
    ungroup()
  ggplot(results, aes(logFC, -log10(`P-value`), colour = direction)) +
    geom_point(aes(alpha = direction != "Not significant")) +
    geom_text_repel(data = labelled_genes, aes(label = symbol), colour = "black") +
    scale_colour_manual(values = c(phenotype_colours, "Not significant" = "grey75"),
                        guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.65, `TRUE` = 1), guide = "none") +
    theme_classic() +
    theme(aspect.ratio = 1,
          axis.line = element_blank(),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.border = element_rect(linewidth = 2, fill = NA),
          panel.grid.major = element_line(linetype = "dotted", colour = "grey85")) +
    labs(x = expression(log[2]~fold~difference), y = expression(-log[10]~p~value),
         title = title)
}

volcano_plots <- imap(bulk_results, ~ plot_volcano(.x, cohort_titles[[.y]]))
ggsave(file.path(figure_directory, "Figure 1A and 1B - Volcano plots.pdf"),
       volcano_plots[["non-COVID"]] | volcano_plots[["COVID"]], width = 12, height = 6)

# Figure 1C: reproducible gene heatmaps ------------------------------------

reproducible_genes <- read_csv(file.path(results_directory, "reproducible_genes.csv"),
                               show_col_types = FALSE)

read_counts <- function(file_name) {
  read_csv(file.path(clean_directory, file_name), show_col_types = FALSE) %>%
    column_to_rownames("ensembl_gene_id") %>%
    as.matrix()
}
cohort_samples <- list(
  "non-COVID" = read_csv(file.path(clean_directory, "noncovid_bulk_subjects.csv"),
                         show_col_types = FALSE),
  "COVID" = read_csv(file.path(clean_directory, "covid_bulk_subjects.csv"),
                     show_col_types = FALSE)
)
cohort_counts <- list(
  "non-COVID" = read_counts("noncovid_bulk_counts.csv"),
  "COVID" = read_counts("covid_bulk_counts.csv")
)

sort_hclust <- function(clustering) as.hclust(dendsort(as.dendrogram(clustering)))

# Expression is DESeq2's variance-stabilized counts, standardized per gene and
# truncated at 2 standard deviations.
draw_reproducible_gene_heatmap <- function(cohort) {
  sample_table <- cohort_samples[[cohort]] %>% column_to_rownames("sample_id")
  stabilized_counts <-
    DESeqDataSetFromMatrix(round(cohort_counts[[cohort]]), colData = sample_table,
                           design = ~ 1) %>%
    vst() %>%
    assay()
  heatmap_matrix <-
    stabilized_counts[reproducible_genes$ensembl_gene_id, ] %>%
    `rownames<-`(reproducible_genes$gene_symbol) %>%
    t() %>% scale() %>% t() %>%
    pmin(2) %>% pmax(-2)
  pheatmap(heatmap_matrix,
           color = mako(99),
           cluster_rows = sort_hclust(hclust(dist(heatmap_matrix))),
           cluster_cols = sort_hclust(hclust(dist(t(heatmap_matrix)))),
           cutree_rows = 2,
           treeheight_row = 10, treeheight_col = 10,
           annotation_col = sample_table %>% dplyr::select(Phenotype = phenotype),
           annotation_colors = list(Phenotype = phenotype_colours),
           show_colnames = FALSE, border_color = NA,
           cellwidth = 10, cellheight = 10,
           main = cohort_titles[[cohort]],
           silent = TRUE)
}

for (cohort in names(cohort_titles)) {
  heatmap <- draw_reproducible_gene_heatmap(cohort)
  pdf(file.path(figure_directory, paste0("Figure 1C - Reproducible genes, ", cohort, ".pdf")),
      width = 14, height = 5)
  grid::grid.draw(heatmap$gtable)
  dev.off()
}

# Reproducible genes in the single-cell RNA-seq -----------------------------

# For each gene, cell type and phenotype: the percentage of cells in which the
# gene is detected (dot size) and its mean expression, log(1 + mean
# normalized expression), standardized per gene across all cell type and
# phenotype groups (colour, on a diverging blue-grey-red scale centred at the
# gene's average; the phenotype colours are not used, so that colour is not
# read as phenotype). Cell types with fewer than 100 cells are left out.
# Genes are grouped by the phenotype in which they are up in bulk RNA-seq.
scseq_seurat <- readRDS(file.path(clean_directory, "scseq_seurat.rds"))

reproducible_gene_symbols <-
  reproducible_genes %>%
  distinct(gene_symbol, direction) %>%
  filter(gene_symbol %in% rownames(scseq_seurat))

# Cell types in order of abundance.
cell_types_shown <-
  scseq_seurat$topcells %>%
  table() %>%
  sort(decreasing = TRUE) %>%
  { names(.)[. >= 100] } %>%
  setdiff("Other")

single_cell_expression <-
  tibble(cell = colnames(scseq_seurat),
         cell_type = scseq_seurat$topcells,
         phenotype = scseq_seurat$phenotype) %>%
  bind_cols(
    LayerData(scseq_seurat, assay = "RNA", layer = "data")[reproducible_gene_symbols$gene_symbol, ] %>%
      as.matrix() %>% t() %>% as_tibble()
  ) %>%
  filter(cell_type %in% cell_types_shown) %>%
  pivot_longer(all_of(reproducible_gene_symbols$gene_symbol),
               names_to = "gene_symbol", values_to = "expression") %>%
  summarize(percent_detected = 100 * mean(expression > 0),
            mean_expression = log1p(mean(expm1(expression))),
            .by = c(gene_symbol, cell_type, phenotype)) %>%
  mutate(scaled_expression = as.numeric(scale(mean_expression)), .by = gene_symbol) %>%
  left_join(reproducible_gene_symbols, by = "gene_symbol") %>%
  mutate(phenotype = factor(phenotype, levels = names(phenotype_labels)),
         cell_type = factor(cell_type, levels = cell_types_shown),
         direction = str_replace(direction, "Up in ", "Up in\n"))

reproducible_gene_dot_plot <-
  ggplot(single_cell_expression,
         aes(phenotype, gene_symbol, size = percent_detected, colour = scaled_expression)) +
  geom_point() +
  facet_grid(direction ~ cell_type, scales = "free_y", space = "free_y",
             labeller = labeller(cell_type = label_wrap_gen(12))) +
  scale_colour_gradient2(name = "Scaled mean\nexpression", low = "#2166AC", mid = "grey85",
                         high = "#B2182B", midpoint = 0,
                         limits = c(-1, 1) * max(abs(single_cell_expression$scaled_expression))) +
  scale_size_area(name = "% cells\ndetected", max_size = 5) +
  scale_x_discrete(labels = c(HypoARDS = "Hypo", HyperARDS = "Hyper")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.border = element_rect(fill = NA, linewidth = 1)) +
  labs(x = "", y = "")

ggsave(file.path(figure_directory, "Reproducible genes by cell type.pdf"),
       reproducible_gene_dot_plot, width = 12.5, height = 5.5)

# UMAPs of the genes up in HyperARDS. The first panel labels the cell types;
# each other panel colours cells by the gene's log-normalized expression
# relative to its 99th percentile among expressing cells (capped at 1), so
# that weakly and strongly expressed genes are equally legible. Expressing
# cells are drawn on top.
umap_cells <-
  Embeddings(scseq_seurat, "umap") %>%
  as_tibble(rownames = "cell") %>%
  dplyr::rename(umap_1 = 2, umap_2 = 3) %>%
  mutate(cell_type = scseq_seurat$topcells[cell])

hyperards_reproducible_genes <-
  reproducible_gene_symbols %>%
  filter(direction == "Up in HyperARDS") %>%
  pull(gene_symbol) %>%
  sort()

umap_expression <-
  umap_cells %>%
  bind_cols(
    LayerData(scseq_seurat, assay = "RNA", layer = "data")[hyperards_reproducible_genes, umap_cells$cell] %>%
      as.matrix() %>% t() %>% as_tibble()
  ) %>%
  pivot_longer(all_of(hyperards_reproducible_genes), names_to = "gene_symbol",
               values_to = "expression") %>%
  mutate(relative_expression = pmin(expression / quantile(expression[expression > 0], 0.99), 1),
         .by = gene_symbol)

umap_theme <- list(
  coord_fixed(),
  theme_void(),
  theme(plot.title = element_text(hjust = 0.5, size = 11))
)

cell_type_umap <-
  ggplot(umap_cells, aes(umap_1, umap_2)) +
  geom_point(colour = "grey80", size = 0.05) +
  geom_label_repel(data = umap_cells %>%
                     summarize(umap_1 = median(umap_1), umap_2 = median(umap_2), .by = cell_type),
                   aes(label = cell_type), size = 2.5, label.size = 0, fill = alpha("white", 0.7),
                   min.segment.length = 0, seed = 1) +
  umap_theme +
  labs(title = "Cell types")

plot_gene_umap <- function(gene) {
  umap_expression %>%
    filter(gene_symbol == gene) %>%
    arrange(relative_expression) %>%
    ggplot(aes(umap_1, umap_2, colour = relative_expression)) +
    geom_point(size = 0.05) +
    scale_colour_viridis_c(name = "Relative\nexpression", limits = c(0, 1)) +
    umap_theme +
    labs(title = gene)
}

ggsave(file.path(figure_directory, "Reproducible genes up in HyperARDS, UMAP.png"),
       wrap_plots(c(list(cell_type_umap), map(hyperards_reproducible_genes, plot_gene_umap)),
                  ncol = 6) +
         plot_layout(guides = "collect"),
       width = 16, height = 9, dpi = 300, bg = "white")

# HyperARDS gene scores. Each cell's score is the mean expression of a set of
# genes up in HyperARDS minus the mean of control genes matched on average
# expression (Seurat::AddModuleScore, 24 expression bins, 100 control genes per
# gene, seed 1); positive scores mean the genes are expressed above matched
# controls. Three sets are scored: all 16 genes, the 7 ribosomal genes
# (cytosolic and mitochondrial ribosomal proteins and RACK1), and the other 9.
#
# Cells within a subject are not independent, so phenotypes are compared at
# the subject level: each point is a subject's mean score, over all of the
# subject's cells or within one cell type. The all-cells score reflects cell
# composition as bulk RNA-seq does; the within-cell-type scores remove it. No
# test is made, with 4 HyperARDS and 5 HypoARDS subjects.
ribosomal_genes <- str_subset(hyperards_reproducible_genes, "^RP[SL]\\d|^MRP[SL]\\d|^RACK1$")
gene_sets <- list(
  "All 16 genes" = hyperards_reproducible_genes,
  "Ribosomal genes" = ribosomal_genes,
  "Other genes" = setdiff(hyperards_reproducible_genes, ribosomal_genes)
)

hyperards_gene_scores <-
  Seurat::AddModuleScore(scseq_seurat, features = unname(gene_sets), assay = "RNA",
                         name = "gene_set_score", seed = 1)[[]] %>%
  as_tibble(rownames = "cell") %>%
  transmute(cell, sample_id, cell_type = topcells,
            phenotype = factor(phenotype, levels = names(phenotype_labels)),
            !!!set_names(syms(paste0("gene_set_score", seq_along(gene_sets))), names(gene_sets))) %>%
  pivot_longer(all_of(names(gene_sets)), names_to = "gene_set", values_to = "score") %>%
  mutate(gene_set = factor(gene_set, levels = names(gene_sets)))

all_gene_scores <-
  hyperards_gene_scores %>%
  filter(gene_set == "All 16 genes") %>%
  inner_join(umap_cells %>% dplyr::select(cell, umap_1, umap_2), by = "cell")

score_limit <- quantile(abs(all_gene_scores$score), 0.99)

score_umap <-
  all_gene_scores %>%
  arrange(abs(score)) %>%
  ggplot(aes(umap_1, umap_2, colour = pmax(pmin(score, score_limit), -score_limit))) +
  geom_point(size = 0.05) +
  facet_wrap(~ phenotype, labeller = as_labeller(phenotype_labels)) +
  scale_colour_gradient2(name = "HyperARDS\ngene score", low = "#2166AC", mid = "grey85",
                         high = "#B2182B", midpoint = 0, limits = c(-1, 1) * score_limit) +
  coord_fixed() +
  theme_void() +
  theme(strip.text = element_text(size = 12, margin = margin(b = 4)))

ggsave(file.path(figure_directory, "HyperARDS gene score by phenotype, UMAP.png"),
       score_umap, width = 11, height = 5, dpi = 300, bg = "white")

# Scores over all of a subject's cells, then within each cell type.
score_groups <- c("All cells", cell_types_shown)
grouped_scores <-
  bind_rows(
    hyperards_gene_scores %>% mutate(group = "All cells"),
    hyperards_gene_scores %>% filter(cell_type %in% cell_types_shown) %>% mutate(group = cell_type)
  ) %>%
  mutate(group = factor(group, levels = score_groups))

subject_mean_scores <-
  grouped_scores %>%
  summarize(score = mean(score), cells = n(), .by = c(sample_id, phenotype, group, gene_set))

score_by_cell_type <-
  ggplot(grouped_scores, aes(phenotype, score)) +
  geom_violin(aes(fill = phenotype), colour = NA, alpha = 0.35, scale = "width") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_point(data = subject_mean_scores, aes(colour = phenotype, size = cells),
             position = position_jitter(width = 0.08, seed = 1)) +
  facet_grid(gene_set ~ group, labeller = labeller(group = label_wrap_gen(12))) +
  scale_fill_manual(values = phenotype_colours, guide = "none") +
  scale_colour_manual(values = phenotype_colours, guide = "none") +
  scale_size_area(name = "Cells from\nsubject", max_size = 4) +
  scale_x_discrete(labels = c(HypoARDS = "Hypo", HyperARDS = "Hyper")) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 1)) +
  labs(x = "", y = "HyperARDS gene score",
       caption = "Violins: all cells. Points: subject means.")

ggsave(file.path(figure_directory, "HyperARDS gene score by cell type and phenotype.pdf"),
       score_by_cell_type, width = 15, height = 10)

# Composition and score together, one point per subject: the percentage of
# the subject's cells that are T or NK cells against the subject's mean
# score over all cells.
subject_composition_scores <-
  hyperards_gene_scores %>%
  summarize(t_nk_percent = 100 * mean(cell_type %in% c("T cell", "NK cell")),
            score = mean(score), cells = n(),
            .by = c(sample_id, phenotype, gene_set))

composition_score_plot <-
  ggplot(subject_composition_scores, aes(t_nk_percent, score, colour = phenotype)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_point(aes(size = cells)) +
  geom_text_repel(aes(label = sample_id), size = 3, show.legend = FALSE, seed = 1) +
  facet_wrap(~ gene_set, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = phenotype_colours, labels = phenotype_labels, name = NULL) +
  scale_size_area(name = "Cells", max_size = 5) +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(fill = NA, linewidth = 1)) +
  labs(x = "T and NK cells (% of subject's cells)",
       y = "Subject mean HyperARDS gene score, all cells")

ggsave(file.path(figure_directory, "HyperARDS gene score and T and NK cell fraction.pdf"),
       composition_score_plot, width = 12, height = 4.5)

# Figures 2A and E2: GSVA scores --------------------------------------------

gsva_path <- file.path(results_directory, "Supplementary Data E3 - GSVA.xlsx")
gsva_tests <- read_csv(file.path(results_directory, "gsva_phenotype_tests.csv"),
                       show_col_types = FALSE)

significance_stars <- function(p_value) {
  case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**", p_value < 0.05 ~ "*", TRUE ~ "ns")
}

plot_gsva_scores <- function(scores, title, test_cohort) {
  stars <- significance_stars(gsva_tests$p_value[gsva_tests$cohort == test_cohort])
  scores %>%
    mutate(phenotype = factor(phenotype, levels = names(phenotype_labels))) %>%
    ggplot(aes(phenotype, `GSVA score`)) +
    geom_boxplot(aes(colour = phenotype), outlier.shape = NA, linewidth = 1) +
    geom_jitter(width = 0.1, size = 3) +
    annotate("segment", x = 1, xend = 2, y = 0.9, yend = 0.9, linewidth = 1) +
    annotate("text", x = 1.5, y = 0.95, label = stars, size = 10) +
    scale_colour_manual(values = phenotype_colours, guide = "none") +
    scale_x_discrete(labels = phenotype_labels) +
    coord_cartesian(ylim = c(-0.75, 1)) +
    theme_classic() +
    theme(aspect.ratio = 1,
          panel.border = element_rect(fill = NA, linewidth = 2),
          axis.line = element_blank(),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 18),
          plot.title = element_text(hjust = 0.5, size = 20)) +
    labs(x = "", y = "GSVA score", title = title)
}

noncovid_gsva <- read_excel(gsva_path, sheet = "S3A - ALI")
covid_gsva <- read_excel(gsva_path, sheet = "S3B - COMET")

ggsave(file.path(figure_directory, "Figure 2A - GSVA scores.pdf"),
       plot_gsva_scores(noncovid_gsva, cohort_titles[["non-COVID"]], "non-COVID") /
         plot_gsva_scores(covid_gsva, cohort_titles[["COVID"]], "COVID, all samples"),
       width = 8, height = 12)

ggsave(file.path(figure_directory, "Figure E2 - GSVA scores, first COVID sample.pdf"),
       covid_gsva %>%
         arrange(subject, vent_day) %>%
         distinct(subject, .keep_all = TRUE) %>%
         plot_gsva_scores("COVID ARDS", "COVID, first sample per subject"),
       width = 12, height = 8)

# Figure 2B: shared bulk pathways ------------------------------------------

# Pathways chosen for display from the pathways enriched in the same direction
# in both bulk cohorts.
pathways_to_plot <- c(
  "NONSENSE MEDIATED DECAY NMD",
  "SELENOAMINO ACID METABOLISM",
  "CELLULAR RESPONSE TO STARVATION",
  "MITOCHONDRIAL TRANSLATION",
  "THE CITRIC ACID TCA CYCLE AND RESPIRATORY ELECTRON TRANSPORT",
  "TCF DEPENDENT SIGNALING IN RESPONSE TO WNT",
  "INTERFERON ALPHA BETA SIGNALING",
  "INTERFERON GAMMA SIGNALING",
  "ESTROGEN DEPENDENT GENE EXPRESSION",
  "DEFECTIVE CFTR CAUSES CYSTIC FIBROSIS",
  "CELLULAR RESPONSE TO HYPOXIA",
  "TCR SIGNALING",
  "DECTIN 1 MEDIATED NONCANONICAL NF KB SIGNALING",
  "CYTOKINE SIGNALING IN IMMUNE SYSTEM",
  "SIGNALING BY THE B CELL RECEPTOR BCR",
  "CELLULAR RESPONSE TO HEAT STRESS",
  "PTEN REGULATION",
  "APOPTOSIS",
  "FC EPSILON RECEPTOR FCERI SIGNALING",
  "NEUTROPHIL DEGRANULATION",
  "GLYCOSAMINOGLYCAN METABOLISM",
  "IMMUNOREGULATORY INTERACTIONS BETWEEN A LYMPHOID AND A NON LYMPHOID CELL",
  "EXTRACELLULAR MATRIX ORGANIZATION",
  "SIGNALING BY GPCR",
  "COLLAGEN BIOSYNTHESIS AND MODIFYING ENZYMES",
  "POTASSIUM CHANNELS"
)

bulk_gsea_shared_pathways <- read_csv(file.path(results_directory, "bulk_gsea_shared_pathways.csv"),
                                      show_col_types = FALSE)
missing_pathways <- setdiff(pathways_to_plot, bulk_gsea_shared_pathways$pathway)
if (length(missing_pathways) > 0) {
  warning("Pathways chosen for Figure 2B but not shared by both cohorts: ",
          paste(missing_pathways, collapse = "; "))
}

shared_pathway_plot <-
  bulk_gsea_shared_pathways %>%
  filter(pathway %in% pathways_to_plot) %>%
  transmute(pathway, ALI = noncovid_nes, COMET = covid_nes) %>%
  pivot_longer(-pathway, names_to = "cohort", values_to = "NES") %>%
  mutate(nes_product = prod(NES) * sign(sum(NES)), .by = pathway) %>%
  mutate(pathway = str_wrap(pathway, indent = 14, width = 43)) %>%
  ggplot(aes(cohort, reorder(pathway, nes_product), size = abs(NES), colour = NES)) +
  geom_point() +
  scale_colour_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(4, 8), guide = "none") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect(linewidth = 2, fill = NA),
        aspect.ratio = 6) +
  labs(x = "", y = "")

ggsave(file.path(figure_directory, "Figure 2B - Shared bulk pathways.pdf"),
       shared_pathway_plot, width = 8, height = 12)

# Figure E3: bulk and single-cell pathways ---------------------------------

# NES is positive for pathways enriched in HyperARDS and negative for those
# enriched in HypoARDS, so it is coloured on a diverging scale from the
# HypoARDS to the HyperARDS colour, through white at zero.
comparison_labels <- c(
  "non-COVID" = "ALI", "COVID" = "COMET", "Alveolar mac" = "AlvMac",
  "Monocyte-derived mac" = "BMDM", "Monocyte" = "Monocytes", "Neutrophil" = "Neutrophils",
  "T cell" = "T cells"
)
sequencing_labels <- c("Bulk RNA-seq" = "Bulk RNASeq", "Single-cell RNA-seq" = "scRNASeq")

bulk_and_single_cell_gsea <- read_csv(file.path(results_directory, "bulk_and_single_cell_gsea.csv"),
                                      show_col_types = FALSE)

bulk_and_single_cell_plot <-
  bulk_and_single_cell_gsea %>%
  mutate(comparison = factor(comparison_labels[comparison], levels = comparison_labels),
         sequencing = sequencing_labels[sequencing],
         pathway = str_wrap(pathway, width = 60)) %>%
  mutate(nes_sum = sum(NES), .by = pathway) %>%
  ggplot(aes(comparison, reorder(pathway, nes_sum), colour = NES, size = abs(NES))) +
  facet_grid(~ sequencing, scales = "free_x", space = "free") +
  geom_point() +
  scale_colour_gradient2(low = phenotype_colours[["HypoARDS"]], mid = "white",
                         high = phenotype_colours[["HyperARDS"]], midpoint = 0) +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.border = element_rect(fill = NA, linewidth = 2),
        strip.background = element_blank()) +
  labs(x = "", y = "")

ggsave(file.path(figure_directory, "Figure E3 - Bulk and single-cell pathways.pdf"),
       bulk_and_single_cell_plot, width = 9,
       height = 2 + 0.4 * n_distinct(bulk_and_single_cell_gsea$pathway))

message("Figures written to ", figure_directory, ":")
message(paste(" ", list.files(figure_directory), collapse = "\n"))
