library(tidyverse)
library(qiime2R)

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
## 3) core feature summaries
################################################################################

## create a .csv file summarizing the taxonomic information available for each of the core sequence features
taxa %>% 
    filter(FeatureID %in% coreFeat_10) %>% 
    select(-Kingdom, -Phylum, -Confidence, -Consensus)
    
  


