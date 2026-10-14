# ARDS Phenotypes - 7 - Figures.R
#
# Draws the transcriptomic figures of the AJRCCM 2025 paper:
#   Figure 1A, 1B  volcano plots of bulk differential expression (Data E1)
#   Figure 1C      heatmaps of the reproducible genes in each bulk cohort
#   Figure 2A      HyperARDS signature GSVA scores by phenotype
#   Figure 2B      Reactome pathways enriched in the same direction in both
#                  bulk cohorts (a curated subset)
#   Figure E2      GSVA scores, first COVID sample per subject
#   Figure E3      the shared bulk pathways that are also enriched in at least
#                  one single-cell type
#   Ligand-receptor figures (script 6 results):
#     the differential interaction network between phenotypes, and a chord
#     diagram of the interactions up in each phenotype (Data E6)
#
# Phenotypes are coloured from the Okabe-Ito palette. Continuous values that
# run from low to high use viridis palettes; Figure E3's NES, centred on zero,
# uses a diverging scale between the phenotype colours. The ligand-receptor
# figures keep the cell type colours of the published figures (CellChat's
# default palette).
# Heatmap dendrograms are ordered by dendsort.
#
# Run from the repository root after scripts 0, 3, 4, 5 and 6. Writes PDFs to
# data/results/ajrccm_2025/figures/.

library(tidyverse)
library(readxl)
library(DESeq2)
library(pheatmap)
library(dendsort)
library(ggrepel)
library(patchwork)
library(viridis)
library(CellChat)

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

# Ligand-receptor figures -------------------------------------------------

cellchat_cell_types <- c("T cell", "Monocyte", "Dendritic cell", "Neutrophil", "Alveolar mac",
                         "Monocyte-derived mac", "NK cell", "B cell")
# The published figures place and colour the cell types in the order of
# cellchat_cell_types, with CellChat's default palette. CellChat itself orders
# them alphabetically, so the figures reorder them.
cell_type_colours <- set_names(scPalette(length(cellchat_cell_types)), cellchat_cell_types)

# Differential interaction network: a circle of cell types whose edges show
# the difference, HyperARDS minus HypoARDS, in the total strength of inferred
# interactions from each cell type to each other.
# Red edges are higher in HyperARDS and blue edges higher in HypoARDS; edge
# width scales with the difference. The merged object from script 6 holds
# HypoARDS first and HyperARDS second, so comparison = c(1, 2).
merged_cellchat <- readRDS(file.path(results_directory, "cellchat_hyperards_hypoards.rds"))

# CellChat 1.5.0 sets a loop angle only on self-loops and leaves the other
# edges missing, which igraph 1.x replaced with 0 but igraph 2.x cannot (its
# default is empty). Restore the default of 0 for this session.
igraph_default_values <- igraph:::i.default.values
igraph_default_values$edge$loop.angle <- 0
assignInNamespace("i.default.values", igraph_default_values, ns = "igraph")
stopifnot(identical(names(merged_cellchat@net), c("HypoARDS", "HyperARDS")),
          setequal(levels(merged_cellchat@idents$joint), cellchat_cell_types))
# netVisual_diffInteraction places cell types in the row order of the
# interaction matrices.
for (phenotype in names(merged_cellchat@net)) {
  merged_cellchat@net[[phenotype]]$weight <-
    merged_cellchat@net[[phenotype]]$weight[cellchat_cell_types, cellchat_cell_types]
}

pdf(file.path(figure_directory, "Ligand-receptor interactions, HyperARDS vs HypoARDS.pdf"),
    width = 6, height = 6)
netVisual_diffInteraction(merged_cellchat, comparison = c(1, 2), measure = "weight",
                          weight.scale = TRUE, color.use = cell_type_colours,
                          title.name = "Interaction strength")
dev.off()

# Chord diagrams of the interactions in Supplementary Data E6 that are up in
# each phenotype, from ligand-expressing to receptor-expressing cell types,
# one page each. netVisual_chord_gene needs only the cell types from its
# CellChat object when the interactions are supplied, so a minimal object
# stands in.
ligand_receptor_path <- file.path(results_directory, "Supplementary Data E6 - Ligand-receptor pairs.xlsx")
cell_type_holder <- new("CellChat", idents = factor(character(0), levels = cellchat_cell_types))
upregulated_sheets <- c(HyperARDS = "6A - Hyperinflammatory ARDS",
                        HypoARDS = "6B - Hypoinflammatory ARDS")

pdf(file.path(figure_directory, "Ligand-receptor pairs up in each phenotype.pdf"),
    width = 10, height = 10)
iwalk(upregulated_sheets, function(sheet, phenotype) {
  netVisual_chord_gene(cell_type_holder, net = as.data.frame(read_excel(ligand_receptor_path, sheet)),
                       color.use = cell_type_colours, lab.cex = 0.5,
                       title.name = paste("Upregulated signaling in", phenotype))
})
dev.off()

message("Figures written to ", figure_directory, ":")
message(paste(" ", list.files(figure_directory), collapse = "\n"))
