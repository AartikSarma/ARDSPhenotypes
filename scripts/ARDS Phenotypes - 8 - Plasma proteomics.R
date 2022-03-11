library(readxl)
library(tidyverse)
library(pheatmap)

colors <- brewer.pal(n = 3, name = "Dark2")
col.hyper <- colors[2]
col.hypo <-  colors[3]
col.control <- colors[1]


ards.subjects <- read_csv("./data/cohort1/cohort1.subjects.csv")
ards.olink <- read_csv("./data/cohort1/cohort1.olink.csv")

olink.qc <- read_csv("./data/cohort1/olink.qc.csv")
qc.included <- 
  olink.qc %>%
  filter(pctmissing < 0.1) %>%
  pull(biomarker)


ards.olink %<>%
  left_join(ards.subjects %>% dplyr::select(sampleid, phenotype)) %>%
  mutate(phenotype = replace_na(phenotype, "Control")) %>%
  filter(`QC Warning` == "Pass") %>%
  dplyr::select(one_of(c("sampleid", "phenotype", qc.included))) 


olink.wilcox <- ards.olink %>%
  pivot_longer(cols = !matches("sampleid") & !matches("phenotype"), names_to = "biomarker") %>%
  dplyr::select(-sampleid) %>%
  group_by(phenotype, biomarker) %>%
  summarize(value = list(value)) %>% 
  pivot_wider(id_cols = "biomarker", names_from = "phenotype") %>%
  group_by(biomarker) %>%
  mutate(p.hypervshypo =wilcox.test(unlist(HyperARDS), unlist(HypoARDS))$p.value,
         p.hypervscontrol =wilcox.test(unlist(HyperARDS), unlist(Control))$p.value,
         p.hypovscontrol =wilcox.test(unlist(HypoARDS), unlist(Control))$p.value) %>%
  ungroup %>%
  mutate(padj.hypervshypo = p.adjust(p.hypervshypo, method = "BH"),
         padj.hypervscontrol =p.adjust(p.hypervscontrol, method = "BH"),
         padj.hypovscontrol =p.adjust(p.hypovscontrol, method = "BH")) 


protein.fdr <- 0.1

olink.wilcox %<>%
  group_by(biomarker) %>%
  mutate(hypervshypo.ratio = median(unlist(HyperARDS))/median(unlist(HypoARDS))) %>%
  mutate(hypervscontrol.ratio = median(unlist(HyperARDS))/median(unlist(Control))) %>%
  mutate(hypovscontrol.ratio = median(unlist(HypoARDS))/median(unlist(Control))) 

olink.wilcox %>%
  group_by(biomarker) %>%
  mutate(hypervshypo = case_when(padj.hypervshypo >= protein.fdr ~ "No diff", 
                                       hypervshypo.ratio > 1 ~ "Hyper",
                                       hypervshypo.ratio < 1 ~ "Hypo"), 
         hypervscontrol= case_when(padj.hypervscontrol >= protein.fdr ~ "No diff", 
                                       hypervscontrol.ratio > 1 ~ "Hyper",
                                       hypervscontrol.ratio < 1 ~ "Control"),
         hypovscontrol = case_when(padj.hypovscontrol >= protein.fdr ~ "No diff", 
                                       hypovscontrol.ratio > 1 ~ "Hypo",
                                       hypovscontrol.ratio < 1 ~ "Control") 
  ) %>% 
  dplyr::select(biomarker, hypervshypo, hypervscontrol, hypovscontrol) %>%
  column_to_rownames("biomarker") -> olink.annot

olink.annot %>% arrange(hypervshypo, desc(hypervscontrol), desc(hypovscontrol)) %>% 
  rownames -> olink.roworder

pdf(height = 30, width = 6, file = "OLink.Wilcox.pdf")
olink.annot %>% rownames_to_column("biomarker") %>%
  pivot_longer(names_to = "comparison", cols = !matches("biomarker"), values_to = "group") %>% 
  mutate(group = factor(group)) %>%
  mutate(comparison = factor(comparison, levels = c("hypervshypo", "hypervscontrol", "hypovscontrol"))) %>%
  group_by(biomarker) %>%
  mutate(y.axisorder = which(biomarker == olink.roworder)) %>%
  ggplot(aes(y = reorder(biomarker, y.axisorder), x = comparison, fill= group)) + 
  geom_tile() + 
  scale_fill_manual(values = c(col.control, col.hyper, col.hypo, col.nodiff)) + 
  theme(aspect.ratio = 72/3, axis.text.x = element_text(angle = 90)) + 
  labs(x = "", y = "")

dev.off()



olink.nm <-
  ards.olink %>% column_to_rownames("sampleid") %>% dplyr::select(-phenotype) %>% scale

olink.breaks <- quantile_breaks(olink.nm)


hyperhypo.nm <- ards.olink %>%
  filter(phenotype %in% c("HyperARDS", "HypoARDS")) %>%
  column_to_rownames("sampleid") %>%
  dplyr::select(-phenotype) %>%
  scale() 


hypercontrol.nm <- ards.olink %>%
  filter(phenotype %in% c("HyperARDS", "Control")) %>%
  column_to_rownames("sampleid") %>%
  dplyr::select(-phenotype) %>%
  scale() 

hypocontrol.nm <- ards.olink %>%
  filter(phenotype %in% c("Control", "HypoARDS")) %>%
  column_to_rownames("sampleid") %>%
  dplyr::select(-phenotype) %>%
  scale() 


col.nodiff <- "grey75"
ann_colors <- 
  list(phenotype = c("Control" = col.control, "HyperARDS" = col.hyper, "HypoARDS" = col.hypo),
       hypervshypo = c("Control" = col.control, "Hyper" = col.hyper, "Hypo" = col.hypo, "No diff" = col.nodiff ),
       hypervscontrol = c("Control" = col.control, "Hyper" = col.hyper, "Hypo" = col.hypo, "No diff" = col.nodiff ),
       hypovscontrol = c("Control" = col.control, "Hyper" = col.hyper, "Hypo" = col.hypo, "No diff" = col.nodiff ))

pdf(file = "./output/Figure 3 - Differences in plasma protein concentrations.pdf", width = 15, height = 20)
library(ggfortify)
pca_res <- prcomp(ards.olink[,3:75], scale = T)
autoplot(pca_res, data = ards.olink, colour = "phenotype") + 
  theme_classic() + 
  geom_point(size = 3, aes(color = phenotype)) + 
  theme(axis.line = element_line(size = 2), 
        axis.title = element_text(size = 20), 
        aspect.ratio = 1) + 
  scale_color_manual(values = c(col.control, col.hyper, col.hypo))


pheatmap(hyperhypo.nm %>% t, rotate = T, 
                breaks = olink.breaks, color = RColorBrewer::brewer.pal(9, "Purples"),
                border_color = NA, 
                clustering_distance_cols = "canberra",
         clustering_distance_rows = "canberra",
         cluster_rows = T,
         show_colnames = F,
         cutree_cols = 2,
         cutree_rows =1,
         annotation_colors = ann_colors,
         uselastdend = F,fontsize = 12,
         annotation_col = olink.phenotypes, 
         annotation_row = olink.annot %>% dplyr::select(hypervshypo), 
         legend = T, annotation_legend = F,cellwidth = 15, cellheight = 15, 
         treeheight_col = 10, treeheight_row = 10, fontsize_row = 10)

olink.phenotypes <- ards.subjects %>% column_to_rownames("sampleid") %>% dplyr::select(phenotype)
pheatmap(hypercontrol.nm %>% t, rotate = T, 
         breaks = olink.breaks, color = RColorBrewer::brewer.pal(9, "Purples"),
         border_color = NA, 
         clustering_distance_cols = "canberra",
         clustering_distance_rows = "canberra",
         cluster_rows = T,
         show_colnames = F,
         cutree_cols = 2,
         cutree_rows =2,
         annotation_colors = ann_colors,
         uselastdend = F,fontsize = 12,
         annotation_col = olink.phenotypes, 
         annotation_row = olink.annot %>% dplyr::select(hypervscontrol), 
         legend = T, annotation_legend = F,cellwidth = 15, cellheight = 15, 
         treeheight_col = 10, treeheight_row = 10, fontsize_row = 10)

pheatmap(hypocontrol.nm %>% t, 
         breaks = olink.breaks, color = RColorBrewer::brewer.pal(9, "Purples"),
         border_color = NA, 
         clustering_distance_cols = "canberra",
         clustering_distance_rows = "canberra",
         cluster_rows = T,
         show_colnames = F,
         cutree_cols = 2,
         cutree_rows =2,
         annotation_colors = ann_colors,
         fontsize = 12,
         annotation_col = olink.phenotypes, 
         annotation_row = olink.annot %>% dplyr::select(hypovscontrol), 
         legend = T, annotation_legend = F,cellwidth = 15, cellheight = 15, 
         treeheight_col = 10, treeheight_row = 10, fontsize_row = 10)




dev.off()

