library(metagenomeSeq)
colors <- brewer.pal(n = 3, name = "Dark2")
col.hyper <- colors[2]
col.hypo <-  colors[3]
col.control <- colors[1]
mngs.ps <- readRDS("./output/ARDS.genus.ps.RDs")
lungData <- mngs.ps %>% 
  prune_samples(sample_sums(.) > 0, .) %>%
  phyloseq_to_metagenomeSeq()


lungData = filterData(lungData, present = 5, depth = 10)
lungData <- wrenchNorm(lungData, condition = lungData$phenotype)
pd <- pData(lungData)
mod <- model.matrix(~0 + phenotype, data = pd)
lungres1 = fitFeatureModel(lungData, mod)


MRcoefs(lungres1,number =  Inf,2) %>% as.data.frame %>%
  arrange(adjPvalues) %>%
  rownames_to_column("genus") %>%
  mutate(genusLabel = case_when(
    adjPvalues < 0.1 ~ genus
  ), 
  logFC = -logFC, 
  logP = -log(adjPvalues, 10), 
  colSig = case_when(
    adjPvalues < 0.1 & logFC > 0 ~ col.hyper, 
    adjPvalues < 0.1 & logFC < 0 ~ col.hypo, 
    TRUE ~ "grey75"
  ), 
  alphaSig = case_when(
    adjPvalues < 0.1 ~ 1, 
    TRUE ~ 0.5
  ))%>% 
  ggplot(aes(x = logFC, y = logP,
             color = colSig, alpha = alphaSig,
             label = genusLabel)) + 
  geom_point() + 
  theme_classic() + 
  scale_color_identity() + 
  scale_alpha_identity() + 
  geom_label_repel() + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 15), 
        axis.line = element_line(size = 1)) + 
  labs(x = bquote("Log_2 fold difference in taxon abundance"), 
       y = expr("-log(p-adjusted)")) -> plot.mngs.diff

plot.mngs.diff
