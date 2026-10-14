# mngs.statistics.r 
# Last updated: November 4, 2021

# Description: This script does this following: 
# 1) Imports data from IDSeq and merges with a sample metadata table
# 2) Uses taxize to import taxnomic data from NCBI and aggregate reads at the taxonomic level of interest
# 3) Applies a negative binomial model to filter out background contamination
# 4) Creates a phyloseq object and calculates alpha/beta diversity metrics

# Required input: sample_taxon_file, sample_overviews_file, mngs_sample_metadata
# Output: 
# 1) mngs.ps: A phyloseq object that contains background-corrected reads mapped
# to a taxon, aggregated at the rank of interest
# 2) mngs.alpha: Alpha diversity metrics for the samples
# 3) mngs.plots.pdf: 
## a) most common taxa: NT_r for taxa vs. sample_erccs for 10 most common taxa
## b) burden of microbial RNA: total NT_r for all included taxa
## c) alpha diversity: violin plots for alpha diversity
## d) beta diversity: dispersion plots for beta diversity estimates

#Load libraries

library(data.table)
library(tidyverse)
library(magrittr)
library(lme4)
library(vegan)
library(phyloseq)
library(taxize)
library(fantaxtic)
library(ggpubr)

ards.subjects <- read_csv("./data/cohort1/cohort1.subjects.csv")
########## Set options here! #######
workingdirectory <- "~/Research/ARDSPhenotypes/"
outputfolder <- "mngs.output"
setwd(workingdirectory)

#Files from IDSeq
sample_taxon_file <- "./data/cohort1/ards_pheno_idseq_counts.csv"
sample_overviews_file <- "./data/cohort1/ards_pheno_idseq_overviews.csv"


# Import IDSeq results
df <- read_csv(sample_taxon_file)
overview <- read_csv(sample_overviews_file)
idseqbatches <- read_csv("./data/cohort1/idseqbatches.csv")

#sample_metadata: A data frame that includes the metadata for all included samples.
# *** SAMPLES THAT ARE NOT IN THIS TABLE WILL BE EXCLUDED *** #
# *** ALL SAMPLES IN sample_data MUST BE BE IN sample_taxon_file *** #
mngs_sample_metadata <- ards.subjects %>%
  filter(IDSeq)

 
dim(mngs_sample_metadata)


#ERCC thresholds for samples to include in background model
background_ercc_max = 1e8
background_ercc_min = 1e5
background_padj_threshold <- 0.05
minsamplesums <- 1

ranklist <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "taxon")
#Set the taxonomic rank to aggregate samples
agg_rank <- ranklist[7] #6 = family, 7 = genus, 8 = taxon

#Microbial kingdoms to include
include_bacteria <- TRUE
include_fungi <- TRUE
include_viruses <- TRUE

####################################


#Clean up the names for downstream analysis
colnames(df) <- make.names(colnames(df))
overview$sampleid %<>% make.names
mngs_sample_metadata$sampleid <- make.names(mngs_sample_metadata$sampleid)


#If every sample in the mngs_sample_metadata object isn't also in the 
#sample_taxon_file, this generates problems in downstream analyses
stopifnot(mngs_sample_metadata$sampleid %in% colnames(df))

## Identify control samples
overview %<>%
  mutate(sample_erccs = total_ercc_reads) %>%
  group_by(sampleid) %>%
  mutate(is_control = str_detect(sampleid, "Control") ) %>%
  mutate(is_control = is_control & sample_erccs > background_ercc_min & sample_erccs < background_ercc_max ) %>%
  ungroup()


overview %<>% 
  dplyr::select(sampleid, sample_name, nonhost_reads, sample_erccs, is_control)

overview %<>% left_join(idseqbatches)

overview %>% filter(is_control) %>% 
  filter(batch == "A" | batch == "B") %>%
  ggplot(aes(x = sample_erccs, y = nonhost_reads, color = batch)) +
  geom_point() + scale_x_continuous(trans = "log10") + 
  scale_y_continuous(trans = "log10")#+  geom_label_repel(aes(label = sample_name))

overview %<>% filter(batch == "B")
#Clean up taxon names
df %<>%
  mutate(taxon =
           str_remove(taxon, "non-species-specific reads in genus ") %>%
           str_remove("\\[") %>%
           str_remove("\\]") %>%
           str_remove("uncultured ")
  ) 

df %<>% dplyr::select(taxon, everything())

df %<>% dplyr::select(one_of(c("taxon", overview$sampleid)))

############################################
# Use this code block to import taxonomy data from NCBI.

# To access the NCBI taxonomy database, you'll need to set up an NCBI API key
# Instructions here: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
# You can then set your ENTREZ_KEY in .Renvrion:
## Run usethis::edit_r_environ()
## Add ENTREZ_KEY = your API key
## Save .Renviron
## Restart R (Shift + Cmd + F10 on a Mac)

#Import taxonomic data from NCBI
## This data import sometimes fails if you try to import too many ranks at once
## The following loops breaks the list of taxa into smaller chunks and tracks
## progress so that you can just restart the loop if the import fails
# taxonomy <- NULL
# counter <-1
# taxonnames <- df$taxon %>% unique
# 
# #If data import error, restart here
# for(i in seq(counter, length(taxonnames), 100)){
#    counter <<- i
#    taxtoget <- taxonnames[i:(i+99)]
#    tmp <- taxize::classification(messages = F,db = "ncbi",sci = taxtoget)
#    taxonomy <<- c(taxonomy, tmp)
# }
# 
# taxonomy<- taxonomy[which(!is.na(taxonomy))]
# 
# save(taxonomy,file =  "mngstaxa.Rda")
#############################################
load(file =  "./data/mngstaxa.Rda")

taxonomy <- rbindlist(taxonomy,idcol = T) %>%
  group_by(.id) %>%
  distinct(rank, .keep_all = T) %>%
  dplyr::select(-id) %>%
  pivot_wider(values_from = "name", names_from = "rank") %>%
  dplyr::rename(taxon = .id) %>%
  dplyr::select(all_of(ranklist)) %>%
  ungroup


# Filter non-microbial reads and
# add taxonomic rank of interest to the sample_taxon_results table

df %<>%
  left_join(taxonomy %>% dplyr::select(taxon, superkingdom, kingdom)) %>%
  filter((superkingdom == "Bacteria" & include_bacteria)|
           (kingdom == "Fungi" & include_fungi) |
           (superkingdom == "Viruses" & include_viruses)) %>%
  dplyr::select(-superkingdom, -kingdom)


df %<>%
  pivot_longer(names_to = "sampleid", cols = !(matches("taxon")),values_to = "NT_r") %>%
  dplyr::select(taxon, matches(agg_rank), everything()) %>%
  group_by(taxon, sampleid) %>%
  summarize(NT_r = sum(NT_r ,na.rm = T))%>%
  left_join(overview) %>%
  dplyr::select(taxon, sampleid, sample_erccs, batch, is_control, NT_r) %>%
  filter(!is.na(taxon)) %>%
  ungroup


#### Background correction model - code from Jack Kamm, adapted from Mick and Kamm, et al. ####
min_reads <- 1
df$adjusted_count <- df$NT_r + min_reads

# # This model takes ~5 hours to complete on an M1 MacBook Pro
bg_model <- glmer.nb(
  formula=adjusted_count ~
    offset(log(sample_erccs)) + (1|taxon),
  data=df[df$is_control,], verbose = 2
)
th <- getME(bg_model, "glmer.nb.theta")

df %>%
  dplyr::mutate(
  #Use sample ERCCs to predict taxon counts
    bg=predict(bg_model, newdata=df, type="response",
               allow.new.levels=TRUE)) %>%
  dplyr::mutate(
  #Calculate p-value for taxon counts in sample vs. background
      pval=pnbinom(adjusted_count, size=th, mu=bg,
                 lower.tail=F)) %>%
  # FWER multiple-testing correction
  dplyr::group_by(sampleid) %>%
  dplyr::mutate(
    padj=p.adjust(pval, method="BH")) %>%
  dplyr::ungroup() ->
  df_with_pvals

gc()
#Remove samples that do not have group metadata
df_with_pvals %<>%
  left_join(mngs_sample_metadata) %>%
  filter(!is.na(phenotype))

saveRDS(df_with_pvals, "./output/mngs.df.batchB.Rds")