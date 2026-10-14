mngs.ps <- readRDS("./output/ARDS.genus.ps.RDs")
colors <- brewer.pal(n = 3, name = "Dark2")
col.hyper <- colors[2]
col.hypo <-  colors[3]
col.control <- colors[1]

set.seed(1)  
#### Alpha diversity ####
tax_table(mngs.ps) %>%
  as.data.frame %>%
  filter(kingdom == "Fungi") %>% 
  rownames() -> fungilist

fungi.ps <- prune_taxa(fungilist, mngs.ps)


tax_table(mngs.ps) %>%
  as.data.frame %>%
  filter(superkingdom == "Bacteria") %>% 
  rownames() -> bacterialist

bacteria.ps <- prune_taxa(bacterialist, mngs.ps)

tax_table(mngs.ps) %>%
  as.data.frame %>%
  filter(
    (superkingdom == "Viruses")
    )%>% 
  rownames() -> viruslist

virus.ps <- prune_taxa(viruslist, mngs.ps)

estimate_richness(mngs.ps, measure = "Shannon") %>% 
  rownames_to_column("sampleid") %>%
  mutate(group = "All taxa")%>%
  bind_rows(
    estimate_richness(bacteria.ps, measure = "Shannon") %>% 
      rownames_to_column("sampleid") %>%
      mutate(group = "Bacteria")
  ) %>%
  bind_rows(
    estimate_richness(fungi.ps, measure = "Shannon") %>% 
      rownames_to_column("sampleid") %>%
      mutate(group = "Fungi")
  ) %>%
  left_join(ards.subjects %>% dplyr::select(sampleid, phenotype)) %>%
  mutate(pheno.col = case_when(
    phenotype == "HyperARDS" ~ col.hyper, 
    phenotype == "HypoARDS" ~ col.hypo
  )) %>%
  ggplot(aes(x = phenotype, y = Shannon, color = pheno.col)) + 
  geom_violin(scale = "width", draw_quantiles = 0.5, size = 1) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1)) + 
  theme_classic() + 
  facet_wrap(~group) + 
  theme(aspect.ratio = 3, axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        axis.line = element_line(size = 1) , 
        strip.text = element_text(size =14),
        strip.background = element_blank()) + 
  scale_color_identity() + 
  stat_compare_means(comparisons = list(c("HyperARDS", "HypoARDS")),) + 
  scale_y_continuous(limits = c(0,3.75)) +
  labs(y = "Shannon diversity index", x = "") -> plot.shannonindex
plot.shannonindex

#### Beta diversity - Exacerbation  ####
exacerbation.samples <- ards.subjects %>%
  filter(IDSeq) %>%
  arrange(sampleid) %>%
  column_to_rownames("sampleid")

ards.subjects %<>%
  group_by(sampleid) %>%
  mutate(pneumonia  = "Pneumonia" %in% c(ali_risk_day1, ali_risk2_day1, ali_risk3_day1)) %>%
  mutate(sepsis  = "Sepsis" %in% c(ali_risk_day1, ali_risk2_day1, ali_risk3_day1)) %>%
  mutate(aspiration  = "Aspiration" %in% c(ali_risk_day1, ali_risk2_day1, ali_risk3_day1)) %>%
  mutate(onpressorssaps = case_when(sampleid == "TA575" ~ 0, TRUE ~ onpressorssaps))

otu_table(mngs.ps) %>% colnames() == rownames(exacerbation.samples)
bcd <- vegdist(otu_table(mngs.ps) %>% t, upper = F, method="bray",na.rm = T)
beta.permanova <-  adonis2(bcd ~ phenotype,
                           data = exacerbation.samples,sqrt.dist = T, by = "terms")

message(paste("The PERMANOVA p-value for Bray-Curtis distance is",beta.permanova$`Pr(>F)`[1] ))
exac.beta.permanova.p <- beta.permanova$`Pr(>F)`[1]
dispersion<-betadisper(bcd, group = exacerbation.samples$phenotype)


dispersion$vectors %>%
  as.data.frame %>%
  dplyr::select(PCoA1, PCoA2) %>%
  rownames_to_column("sampleid") %>%
  left_join(ards.subjects) %>%
  group_by(phenotype) %>% 
  mutate(pc1end = mean(PCoA1), pc2end = mean(PCoA2)) %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = phenotype)) + 
  geom_point(size = 3) + 
  geom_segment(aes(xend = pc1end, yend = pc2end), size = 0.1) +
  stat_conf_ellipse(size = 0.5) + 
  geom_label(aes(label = phenotype, x = pc1end, y = pc2end)) +
  theme_classic() +
  scale_color_manual(values = c(col.hyper, col.hypo)) + 
  labs(caption = paste("PERMANOVA p-value =", exac.beta.permanova.p)) +
  theme(aspect.ratio = 1, axis.line = element_line(size = 1), legend.position = "none") -> plot.beta.alltaxa

plot.beta.alltaxa

## Bacterial

bacteria.ps %<>% prune_samples(sample_sums(.) > 0, .)
bcd <- vegdist(otu_table(bacteria.ps) %>% t, 
               upper = F,
               method="bray",
               na.rm = T)
sample_names(bacteria.ps)
bacteria.samples <- 
  ards.subjects %>%
  filter(sampleid %in% sample_names(bacteria.ps))
beta.permanova <-  adonis2(bcd ~ phenotype,
                           data = bacteria.samples, sqrt.dist = T)

message(paste("The PERMANOVA p-value for Bray-Curtis distance is",beta.permanova$`Pr(>F)`[1] ))
exac.beta.permanova.p <- beta.permanova$`Pr(>F)`[1]
dispersion<-betadisper(bcd, group = bacteria.samples$phenotype)

dispersion$vectors %>%
  as.data.frame %>%
  dplyr::select(PCoA1, PCoA2) %>%
  rownames_to_column("sampleid") %>%
  left_join(ards.subjects) %>%
  group_by(phenotype) %>% 
  mutate(pc1end = mean(PCoA1), pc2end = mean(PCoA2)) %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = phenotype)) + 
  geom_point(size = 3) + 
  geom_segment(aes(xend = pc1end, yend = pc2end), size = 0.1) +
  stat_conf_ellipse(size = 0.5) + 
  geom_label(aes(label = phenotype, x = pc1end, y = pc2end)) +
  theme_classic() +
  scale_color_manual(values = c(col.hyper, col.hypo)) + 
  labs(caption = paste("PERMANOVA p-value =", exac.beta.permanova.p)) +
  theme(aspect.ratio = 1, axis.line = element_line(size = 1), legend.position = "none") -> plot.beta.bacteria

## Fungi

fungi.ps %<>% prune_samples(sample_sums(.) > 0, .)
bcd <- vegdist(otu_table(fungi.ps) %>% t, 
               upper = F,
               method="bray",
               na.rm = T)
sample_names(fungi.ps)
fungi.samples <- 
  ards.subjects %>%
  filter(sampleid %in% sample_names(fungi.ps))
beta.permanova <-  adonis2(bcd ~ phenotype,
                           data = fungi.samples, sqrt.dist = T)

message(paste("The PERMANOVA p-value for Bray-Curtis distance is",beta.permanova$`Pr(>F)`[1] ))
exac.beta.permanova.p <- beta.permanova$`Pr(>F)`[1]
dispersion<-betadisper(bcd, group = fungi.samples$phenotype)

dispersion$vectors %>%
  as.data.frame %>%
  dplyr::select(PCoA1, PCoA2) %>%
  rownames_to_column("sampleid") %>%
  left_join(ards.subjects) %>%
  group_by(phenotype) %>% 
  mutate(pc1end = mean(PCoA1), pc2end = mean(PCoA2)) %>%
  ggplot(aes(x = PCoA1, y = PCoA2, color = phenotype)) + 
  geom_point(size = 3) + 
  geom_segment(aes(xend = pc1end, yend = pc2end), size = 0.1) +
  stat_conf_ellipse(size = 0.5) + 
  geom_label(aes(label = phenotype, x = pc1end, y = pc2end)) +
  theme_classic() +
  scale_color_manual(values = c(col.hyper, col.hypo)) + 
  labs(caption = paste("PERMANOVA p-value =", exac.beta.permanova.p)) +
  theme(aspect.ratio = 1, axis.line = element_line(size = 1), legend.position = "none") -> plot.beta.fungi