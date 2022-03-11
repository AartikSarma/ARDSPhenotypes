setwd("~/Research/ARDSPhenotypes/")

col.palette <- RColorBrewer::brewer.pal(3, "Dark2")
col.hyper <- col.palette[2]
col.hypo <- col.palette[3]
col.control <- col.palette[1]
col.nodiff <- "grey75"

ards.subjects <- read_csv("./data/cohort1/cohort1.subjects.csv")
ards.rna <- read_csv("./data/cohort1/cohort1.rnaseq.csv")


ensembltohgnc <- ards.rna %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol)


ards.subjects %<>%
  mutate(scaleage = scale(agephi)) %>%
  filter(RNASeq)

ards.rna %<>%
  column_to_rownames("ensembl_gene_id") %>%
  dplyr::select(one_of(ards.subjects$sampleid))


#Exp 0: Compare all ARDS subjects to controls
exp0.subjects <- 
  ards.subjects %>%
  filter(RNASeq) %>%
  mutate(pheno2 = phenotype, 
         phenotype = case_when(phenotype == "Control" ~ "Control", 
                               TRUE ~ "ARDS")) %>%
  mutate(phenotype = factor(phenotype, levels = c("Control", "ARDS")))

exp0.rna <- 
  ards.rna %>%
  dplyr::select(one_of(exp0.subjects$sampleid))

exp0.dds <- 
  DESeqDataSetFromMatrix(
    colData = exp0.subjects, 
    countData = exp0.rna, 
    design = ~RNABatch + gender1 + scaleage + phenotype
  )

exp0.vst <- vst(exp0.dds)
exp0.pca <- plotPCA(exp0.vst, intgroup = "pheno2")

#Exp 1: Compare gene expression between phenotypes
exp1.subjects <- 
  ards.subjects %>%
  filter(RNASeq, phenotype != "Control") %>%
  mutate(phenotype = factor(phenotype, levels = c("HypoARDS", "HyperARDS")))

exp1.rna <- 
  ards.rna %>%
  dplyr::select(one_of(exp1.subjects$sampleid))

exp1.dds <- 
  DESeqDataSetFromMatrix(
    colData = exp1.subjects, 
    countData = exp1.rna, 
    design = ~RNABatch + gender1 + scaleage + phenotype
  )

exp1.vst <- vst(exp1.dds)

#Exp 2: Compare gene expression between hyperinflammatory ARDS and controls
exp2.subjects <- 
  ards.subjects %>%
  filter(RNASeq, phenotype != "HypoARDS") %>%
  mutate(phenotype = factor(phenotype, levels = c("Control", "HyperARDS")))

exp2.rna <- 
  ards.rna %>%
  dplyr::select(one_of(exp2.subjects$sampleid))

exp2.dds <- 
  DESeqDataSetFromMatrix(
    colData = exp2.subjects, 
    countData = exp2.rna, 
    design = ~RNABatch + gender1 + scaleage + phenotype
  )

#Exp 3: Compare gene expression between hypoinflammatory ARDS and controls
exp3.subjects <- 
  ards.subjects %>%
  filter(RNASeq, phenotype != "HyperARDS") %>%
  mutate(phenotype = factor(phenotype, levels = c("Control", "HypoARDS")))

exp3.rna <- 
  ards.rna %>%
  dplyr::select(one_of(exp3.subjects$sampleid))

exp3.dds <- 
  DESeqDataSetFromMatrix(
    colData = exp3.subjects, 
    countData = exp3.rna, 
    design = ~RNABatch + gender1 + scaleage + phenotype
  )

# Run differential expression analyses
exp0.dds %<>% DESeq(parallel = T)
exp1.dds %<>% DESeq(parallel = T)
exp2.dds %<>% DESeq(parallel = T)
exp3.dds %<>% DESeq(parallel = T)

#deseqresults is a convenience function in rnaseqscripts.r
#The function takes a DESeq object as input
#and returns the following in a data frame: 
#1) Shrunken LFC estimates 
#2) HGNC symbols
#3) Independent-hypothesis weighted FDR
#By default, the function calculates these values for 
#all available contrasts in the last term of the design formula
#and uses apelgm with adaptive prior for FC estimates
#In addition, the function saves the results in a new folder in the working directory
setwd("output")
dir.create("deseqresults")
setwd("deseqresults")
deseqresults(exp0.dds)
deseqresults(exp1.dds)
deseqresults(exp2.dds)
deseqresults(exp3.dds)
setwd("..")
setwd("..")

#How many genes were differentially expressed between groups? 
exp0.ndeg <- exp0.ARDSvsControl.ape.adapt %>%
  filter(fdr < 0.1, abs(log2FoldChange) > 0.5) %>% nrow
exp1.ndeg <- exp1.HyperARDSvsHypoARDS.ape.adapt %>%
  filter(fdr < 0.1, abs(log2FoldChange) > 0.5) %>% nrow
exp2.ndeg <- exp2.HyperARDSvsControl.ape.adapt %>%
  filter(fdr < 0.1, abs(log2FoldChange) > 0.5) %>% nrow
exp3.ndeg <- exp3.HypoARDSvsControl.ape.adapt %>%
  filter(fdr < 0.1, abs(log2FoldChange) > 0.5) %>% nrow

message(exp0.ndeg, " genes (",round(exp0.ndeg/nrow(exp0.rna) * 100 ), "% of all protein-coding genes) were differentially expressed between ARDS subjects and controls.")
message(exp1.ndeg, " genes (",round(exp1.ndeg/nrow(exp1.rna) * 100 ), "% of all protein-coding genes) were differentially expressed between ARDS phenotypes.")
message(exp2.ndeg, " genes (",round(exp2.ndeg/nrow(exp2.rna) * 100 ), "% of all protein-coding genes) were differentially expressed between hyperinflammatory ARDS and controls.")
message(exp3.ndeg, " genes (",round(exp3.ndeg/nrow(exp3.rna) * 100 ), "% of all protein-coding genes) were differentially expressed between hypoinflammatory ARDS and controls.")

#Write results to working directory to export to IPA
resultsforipa <- ls(pattern = ".ape.adapt")
setwd("output")
dir.create("exporttoipa")
setwd("exporttoipa")
for (result in resultsforipa){
  filename <- result %>% str_remove(".ape.adapt") %>% paste(".csv", sep = "")
  get(result) %>% dplyr::select(ensembl_gene_id, log2FoldChange, pvalue, fdr) %>%  write_csv(file = filename)
}
setwd("..")
setwd("..")

library(GSVA)
library(sva)
#GSVA compared to experimental signatures of gene expression
gsva.ardssigs <- gmtPathways("./data/geosignatures/ardssigs.gmt")

ards.combat <- #Batch correct count data before GSVA
  ComBat_seq(as.matrix(ards.rna), ards.subjects %>% filter(RNASeq) %>% pull(RNABatch))

#Calculate GSVA scores for each set of genes upregulated in experimental models
ards.gsva <- ards.combat %>%
  as.data.frame %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(ensembltohgnc) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  distinct(hgnc_symbol, .keep_all = T) %>%
  column_to_rownames("hgnc_symbol") %>%
  dplyr::select(-ensembl_gene_id) %>%
  as.matrix %>% gsva(gset.idx.list = gsva.ardssigs,
                     method = "gsva",
                     kcdf = "Poisson" #Important! KCDF must be set to Poisson for count data
                     ) 

#Calculate contrasts for GSVA scores using standard limma pipeline
library(limma)
mod <- model.matrix(~ factor(ards.subjects %>%
                               filter(RNASeq) %>%
                               pull(phenotype)))
colnames(mod) <- c("Control", "HyperARDS", "HypoARDS")

fit <- lmFit(ards.gsva, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.1)

contrast.matrix <- makeContrasts(HyperARDS-Control, HypoARDS-Control, levels=mod)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE) 
res2 <- decideTests(fit2, p.value=0.1)
summary(res2)

tt.hyper <- topTable(fit2, coef = 1, number = Inf) %>%  rownames_to_column("gseModel") %>% mutate(pheno.col = col.hyper)
tt.hypo <- topTable(fit2, coef = 2, number = Inf) %>%   rownames_to_column("gseModel") %>%mutate(pheno.col = col.hypo)

tt.hyper %>%
  bind_rows(tt.hypo) %>% 
  mutate(gseModel = gseModel %>%
           str_remove("_UP") %>%
           str_replace_all("_", " ")) %>%
  mutate(is_sig = adj.P.Val < 0.1) %>%
  group_by(gseModel) %>%
  mutate(sumFC = sum(logFC)) %>%
  ggplot(aes(x = logFC,
             y = reorder(gseModel, sumFC))) +
  geom_vline(xintercept = 0, color = col.control) +
  geom_point(aes(color = pheno.col, 
                 alpha = is_sig, 
                 shape = is_sig), size = 4) + 
  scale_color_identity() + 
  scale_shape_manual(values = c(16,18)) +
  scale_alpha_manual(values = c(0.35, 1))+ 
  ggforestplot::geom_stripes(aes()) +
  theme_classic() + 
  scale_x_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75), 
                     labels=  c(-0.25, 0, 0.25, 0.5, 0.75)) + 
  labs(x = "GSVA score vs. control", 
       y = "") + 
  theme(axis.line = element_line(size = 1), 
        legend.position = "none", 
        aspect.ratio =1.61) -> plot.gsvascores
