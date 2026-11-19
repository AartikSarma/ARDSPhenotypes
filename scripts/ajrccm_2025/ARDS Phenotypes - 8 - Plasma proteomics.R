# ARDS Phenotypes - 8 - Plasma proteomics.R
#
# Compares plasma protein concentrations between HyperARDS and HypoARDS in the
# ALI (non-COVID) cohort (Figure 4, Supplementary Data E7). Ported from
# scripts/preprint_2022/ARDS Phenotypes - 8 - Plasma proteomics.R, keeping its
# quality control and test but dropping the control subjects, which the
# published paper does not include.
#
# Plasma was assayed on the Olink 96-protein panel, which reports normalized
# protein expression (NPX) on a log2 scale. As in the paper's Methods:
#   - samples carrying an Olink quality control warning are excluded;
#   - proteins measurable in fewer than 90% of samples are excluded, using the
#     missingness recorded in olink.qc.csv;
#   - NPX is compared between phenotypes by a Wilcoxon rank-sum test, and
#     proteins with a Benjamini-Hochberg FDR < 0.1, adjusted over the proteins
#     that pass quality control, differ between phenotypes.
# The effect size is the difference in mean NPX, HyperARDS minus HypoARDS,
# which is a log2 fold difference because NPX is already log2.
#
# The paper reports 25 subjects (5 HyperARDS, 20 HypoARDS) and 73 proteins
# after quality control, 28 of them higher in HyperARDS; the script stops if
# the subject or protein counts differ. (The Methods section says 76 proteins;
# the Results, Figure 4 and the preprint code all give 73.)
#
# The Olink data are not public but are available on
# reasonable request, so run_all.R does not run this script. It expects three
# files in the git-ignored folder data/olink/:
#   cohort1.subjects.csv  sampleid, phenotype (HyperARDS, HypoARDS; controls
#                         have no phenotype)
#   cohort1.olink.csv     one row per sample: sampleid, `QC Warning`, and one
#                         NPX column per protein
#   olink.qc.csv          biomarker, pctmissing (fraction of samples in which
#                         the protein could not be measured)
#
# Run from the repository root. Writes to data/results/ajrccm_2025/:
#   Supplementary Data E7 - OLink Proteomics Results.xlsx
#   figures/Figure 4 - Plasma proteomics volcano plot.pdf

library(tidyverse)
library(ggrepel)

olink_directory <- file.path(here::here(), "data", "olink")
results_directory <- file.path(here::here(), "data", "results", "ajrccm_2025")
figure_directory <- file.path(results_directory, "figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

fdr_threshold <- 0.1
maximum_fraction_missing <- 0.1
# Okabe-Ito vermillion and blue, as in script 7.
phenotype_colours <- c(HyperARDS = "#D55E00", HypoARDS = "#0072B2")

input_paths <- file.path(olink_directory,
                         c("cohort1.subjects.csv", "cohort1.olink.csv", "olink.qc.csv"))
if (!all(file.exists(input_paths))) {
  stop("The Olink data are not public; this script needs ",
       paste(basename(input_paths[!file.exists(input_paths)]), collapse = ", "),
       " in ", olink_directory, ".")
}

# Data ----------------------------------------------------------------------

ards_subjects <- read_csv(input_paths[1], show_col_types = FALSE) %>%
  filter(phenotype %in% names(phenotype_colours)) %>%
  distinct(sampleid, phenotype)
olink_measurements <- read_csv(input_paths[2], show_col_types = FALSE)
olink_quality_control <- read_csv(input_paths[3], show_col_types = FALSE)

proteins_passing_qc <- olink_quality_control %>%
  filter(pctmissing < maximum_fraction_missing) %>%
  pull(biomarker)

# Inner join: control samples, which have no ARDS phenotype, drop out here.
protein_expression <- olink_measurements %>%
  filter(`QC Warning` == "Pass") %>%
  inner_join(ards_subjects, by = "sampleid") %>%
  select(sampleid, phenotype, all_of(proteins_passing_qc)) %>%
  pivot_longer(all_of(proteins_passing_qc), names_to = "biomarker", values_to = "npx")

subjects_per_phenotype <- protein_expression %>%
  distinct(sampleid, phenotype) %>%
  count(phenotype) %>%
  deframe()
message("Subjects passing Olink QC: ",
        paste(names(subjects_per_phenotype), subjects_per_phenotype, sep = " = ", collapse = ", "))
message("Proteins passing QC: ", length(proteins_passing_qc))
stopifnot(identical(subjects_per_phenotype[c("HyperARDS", "HypoARDS")],
                    c(HyperARDS = 5L, HypoARDS = 20L)),
          length(proteins_passing_qc) == 73)

# Wilcoxon tests --------------------------------------------------------------

# wilcox.test drops missing NPX values; n_ columns count the values tested.
proteomics_results <- protein_expression %>%
  summarize(
    n_hyperards = sum(!is.na(npx[phenotype == "HyperARDS"])),
    n_hypoards = sum(!is.na(npx[phenotype == "HypoARDS"])),
    mean_npx_hyperards = mean(npx[phenotype == "HyperARDS"], na.rm = TRUE),
    mean_npx_hypoards = mean(npx[phenotype == "HypoARDS"], na.rm = TRUE),
    p_value = wilcox.test(npx[phenotype == "HyperARDS"], npx[phenotype == "HypoARDS"])$p.value,
    .by = biomarker
  ) %>%
  mutate(log2_fold_difference = mean_npx_hyperards - mean_npx_hypoards,
         adjusted_p_value = p.adjust(p_value, method = "BH"),
         higher_in = case_when(adjusted_p_value >= fdr_threshold ~ "Not significant",
                               log2_fold_difference > 0 ~ "HyperARDS",
                               TRUE ~ "HypoARDS")) %>%
  arrange(p_value)

writexl::write_xlsx(
  list("E7 - OLink Proteomics" = proteomics_results),
  file.path(results_directory, "Supplementary Data E7 - OLink Proteomics Results.xlsx"))

message("Proteins differing between phenotypes (FDR < ", fdr_threshold, "):")
print(count(proteomics_results, higher_in))

# Figure 4: volcano plot ------------------------------------------------------

# Every protein that differs between phenotypes is labelled.
proteomics_volcano <-
  ggplot(proteomics_results,
         aes(log2_fold_difference, -log10(adjusted_p_value), colour = higher_in)) +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.5) +
  geom_text_repel(data = filter(proteomics_results, higher_in != "Not significant"),
                  aes(label = biomarker), colour = "black", size = 3.5, max.overlaps = Inf) +
  scale_colour_manual(values = c(phenotype_colours, "Not significant" = "grey75"),
                      guide = "none") +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        panel.border = element_rect(linewidth = 2, fill = NA),
        panel.grid.major = element_line(linetype = "dotted", colour = "grey85")) +
  labs(x = expression(log[2]~fold~difference~"in"~mean~NPX),
       y = expression(-log[10]~adjusted~p~value))

ggsave(file.path(figure_directory, "Figure 4 - Plasma proteomics volcano plot.pdf"),
       proteomics_volcano, width = 7, height = 7)
