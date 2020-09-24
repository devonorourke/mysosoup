library(tidyverse)
library(qiime2R)
library(scales)

## function for plot theme:
theme_devon <- function () {
  theme_bw(base_size=12, base_family="Avenir") %+replace%
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA),
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}


################################################################################
## 1a) import read counts and taxa information
## 1b) import core-features featureIDs (see: https://github.com/devonorourke/mysosoup/blob/master/docs/diversity_workflow.md)
################################################################################
## add taxonomy information
## amend import of taxa info
taxa <- read_csv(file="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/taxonomy/filtd_tax_dataframe_ALL.csv")
colnames(taxa)[1] <- "FeatureID"

## add read count data
## download from github:
#download.file(url = 'https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/Mangan.clust_p985_table_Filtd_min10k.qza')
#coretable <- read_qza('Mangan.clust_p985_table_Filtd_min10k.qza')
## or run from local:
coretable <- read_qza('~/github/mysosoup/data/qiime_qza/Mangan.clust_p985_table_Filtd_min10k.qza')
coredata <- as.data.frame(coretable$data) %>% 
  mutate(FeatureID = row.names(.)) %>% 
  pivot_longer(-FeatureID, names_to = "SampleID", values_to = "Reads") %>% 
  filter(Reads > 0)
colnames(coredata)[1] <- "FeatureID"
rm(coretable)

## import core featureIDs
corefeatIDs <- read_csv(file="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/text_tables/core_features/featureIDs_core10.txt")
colnames(corefeatIDs)[1] <- "FeatureID"

## filter to retain only those core-10% featureIDs
feature_taxa <- taxa %>% filter(FeatureID %in% corefeatIDs$FeatureID)


################################################################################
## 2) data explroation justifying using the 10% threshold for core features
################################################################################
## how many unique samples do we have?
coredata %>% distinct(SampleID) %>% nrow()
## 196 samples

## find FeatureIDs in at least...
## 5% of samples (>= 10 samples)
coreFeat_05 <- coredata %>% 
  group_by(FeatureID) %>% tally() %>%
  filter(n >= 10) %>% select(FeatureID) %>% pull()
  ## 120 features in 5% of samples

## 10% of samples (>= 20 samples)
coreFeat_10 <- coredata %>% 
  group_by(FeatureID) %>% tally() %>%
  filter(n >= 20) %>% select(FeatureID) %>% pull()
  ## 54 features in 10% of samples

## at least two or more samples:
coreFeat_noSingleton <- coredata %>% 
  group_by(FeatureID) %>% tally() %>%
  filter(n >1) %>% select(FeatureID) %>% pull()
  ## 523 features in sequence features present in at least 2 samples

## how many distinct taxonomic Orders? Families? OTUs (i.e. FeatureIDs)?
  ## across entire dataset?
  taxa %>% distinct(Order) %>% nrow()
  taxa %>% distinct(Order, Family) %>% nrow()
  taxa %>% distinct(FeatureID) %>% nrow()
  ## 19 Orders, 175 Families, 1124 OTUs detected in entire dataset

  ## in at least 2 or more samples?
  taxa %>% filter(FeatureID %in% coreFeat_noSingleton) %>% distinct(Order) %>% nrow()
  taxa %>% filter(FeatureID %in% coreFeat_noSingleton) %>% distinct(Order, Family) %>% nrow()
  taxa %>% filter(FeatureID %in% coreFeat_noSingleton) %>% distinct(FeatureID) %>% nrow()
  ## 17 Orders, 107 Families, 523 OTUs detected in at least 2 or more samples
  
  ## across just those OTUs in at least 5% of samples?
  taxa %>% filter(FeatureID %in% coreFeat_05) %>% distinct(Order) %>% nrow()
  taxa %>% filter(FeatureID %in% coreFeat_05) %>% distinct(Order, Family) %>% nrow()
  taxa %>% filter(FeatureID %in% coreFeat_05) %>% distinct(FeatureID) %>% nrow()
  ## 10 Orders, 35 Families, 120 OTUs
  
  ## across just those OTUs in at least 10% of samples?
  taxa %>% filter(FeatureID %in% coreFeat_10) %>% distinct(Order) %>% nrow()
  taxa %>% filter(FeatureID %in% coreFeat_10) %>% distinct(Order, Family) %>% nrow()
  taxa %>% filter(FeatureID %in% coreFeat_10) %>% distinct(FeatureID) %>% nrow()
  ## 8 Orders, 21 Families, 54 OTUs


## sanity check: 
  ## for those 5% and/or 10% core features, how many samples have these OTUs?
  ## plot shows number of samples for each OTU among these core features
  ## red and blue dotted line shows the 5% and 10% thresholds, respectively
coredata %>% 
  filter(FeatureID %in% coreFeat_05) %>% 
  group_by(FeatureID) %>% 
  tally() %>% 
  mutate(FeatureID = fct_reorder(FeatureID, n)) %>% 
  ggplot(aes(FeatureID, n)) +
  geom_col() +
  geom_hline(yintercept = 10, linetype="dotted", color="red") +
  geom_hline(yintercept = 20, linetype="dotted", color="blue")
  ## sticking with 10% value - moost of these are present in more than 30-50 samples
  ## whereas the 5% threshold, you have more than half of the samples with less than 20 samples
  ## ... that could be driven by some samples only collected on a certain month/site, so they're not really "core"
  
## how many of these are in each Order? Family?
  ## across entire dataset?
  taxa %>% group_by(Order) %>% tally()
  ## in at least 2 or more samples?
  taxa %>% filter(FeatureID %in% coreFeat_noSingleton) %>% group_by(Order) %>% tally()
  ## across just those OTUs in top 5% of samples?
  taxa %>% filter(FeatureID %in% coreFeat_05) %>% group_by(Order) %>% tally()
  ## across just those OTUs in top 10% of samples?
  feature_taxa %>% group_by(Order) %>% tally()
  


################################################################################
## 3) alldata vs. core feature summaries
################################################################################
## merge data for calculations
  coredataNtaxa <- merge(coredata, taxa, by="FeatureID", all.x=TRUE)
  ## note we lose 36 OTUs between the coredata and the full taxa... 
  ## this is because we filtered out a few samples with less than 10k reads, and those OTUs were derived from those samples
  ## these samples aren't important to our analyses because they are only present in a handful of low-sequenced samples
  
## Table 1 for paper:
lengthSamples = n_distinct(coredataNtaxa$SampleID)
lengthTaxa = n_distinct(coredataNtaxa$FeatureID)
coreOrder <- coredataNtaxa %>% filter(FeatureID %in% coreFeat_10) %>% distinct(Order) %>% pull()

coreOrderSumry <- coredataNtaxa %>% 
  filter(Order %in% coreOrder) %>%
  group_by(Order) %>%
  summarise(nSamples = n_distinct(SampleID),
            nTaxa = n_distinct(FeatureID)) %>% 
  mutate(pSamples = nSamples / lengthSamples,
         pTaxa = (nTaxa / lengthTaxa)) %>% 
  mutate(pTaxa = round(pTaxa, 3) * 100,
         pSamples = round(pSamples, 3) * 100) %>% 
  select(-nSamples, -nTaxa)
write_csv(coreOrderSumry, path="~/github/mysosoup/data/taxonomy/coreOrderSumry.csv")  
  
## Table 2 for paper:
coreOTUsumry <- coredataNtaxa %>% 
  filter(FeatureID %in% coreFeat_10) %>% 
  select(-Kingdom, -Phylum, -Confidence, -Consensus) %>% 
  group_by(Class, Order, Family, Genus, Species) %>%
  summarise(nSamples = n_distinct(SampleID)) %>% 
  arrange(Order, -nSamples)
write_csv(coreOTUsumry, path="~/github/mysosoup/data/taxonomy/coreOTU_summary.csv")

  
## Other exploratory data to retain:  
##  First, create a simple summary file showing the featureIDs taxonomic information
### create a .csv file summarizing the taxonomic information available for each of the core sequence features
coreTaxaSumry <- taxa %>% 
  filter(FeatureID %in% coreFeat_10) %>% 
  select(-Kingdom, -Phylum, -Confidence, -Consensus) %>% 
  arrange(Class, Order, Family, Genus, Species)
write_csv(coreTaxaSumry, path="~/github/mysosoup/data/taxonomy/coreTaxa_summary.csv")    

allTaxaSumry <- taxa %>% 
  select(-Kingdom, -Phylum, -Confidence, -Consensus) %>% 
  arrange(Class, Order, Family, Genus, Species)
write_csv(allTaxaSumry, path="~/github/mysosoup/data/taxonomy/allTaxa_summary.csv")

## Summarize total number of reads in dataset, and total number of samples in dataset 
sumAllReads <- sum(coredataNtaxa$Reads)
nAllSamples <- n_distinct(coredataNtaxa$SampleID)

perSampSumry <- coredataNtaxa %>% 
  group_by(FeatureID, Class, Order, Family, Genus, Species, Classifier) %>% 
  summarize(nDetections = n(),
            nReads = sum(Reads)) %>% 
  mutate(pDetections = nDetections / nAllSamples,
         pReads = nReads / sumAllReads) %>% 
  arrange(Order, -pDetections)
write_csv(perSampSumry, path="~/github/mysosoup/data/taxonomy/perSample_Detections_SeqAbundance_summary.csv")

## how many OTUs are present in just a single sample?
perSampSumry %>% 
  filter(nDetections <= 2) %>% 
  nrow()
  ungroup() %>% 
  summarise(meanReads = mean(nReads))



#################################
## unused code
#################################
# CoreOrder <- taxa %>% filter(FeatureID %in% coreFeat_10) %>% distinct(Order) %>% pull()
# 
# ggplot() +
#   geom_point(data = perSampSumry %>% filter(!FeatureID %in% coreFeat_10), 
#              aes(x=nDetections, y=pReads), color="gray50") +
#   geom_point(data = perSampSumry %>% filter(FeatureID %in% coreFeat_10), 
#              aes(x=nDetections, y=pReads, color=Order)) +
#   scale_x_continuous(trans="log2", labels = comma) +
#   scale_y_continuous(trans="log2", labels = comma) +
#   theme_devon() +
#   labs(x="fraction of samples with sequence variant detected", 
#        y="fraction of all sequences assigned to sequence variant",
#        color="arthropod\norder")
# 
# #ggplot(data=perSampSumry %>% filter(FeatureID %in% coreFeat_10), 
# ggplot(data=perSampSumry, 
#        aes(x=Order, y=pDetections)) +
#   geom_jitter(width = 0.2) +
#   labs(y="fraction of samples detected", x="") +
#   scale_y_continuous(trans = "log2", breaks = c(0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64)) +
#   theme_devon() +
#   theme(axis.text.x = element_text(angle=22.5, hjust=1))
# #facet_grid(~ Order, space = "free_x", shrink = TRUE, scales = "free_x")
