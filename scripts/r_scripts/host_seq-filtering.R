library(tidyverse)
library(reshape2)
library(qiime2R)
library(scales)

## import host .qza file:
# notrun: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan_BatsOnly_ASVtable.qza", "host.qza")
## convert .qza to matrix, then convert wide-format matrix to long-format data.frame object
tableimport.function <- function(table){
  featuretable <- read_qza(table)
  mat.tmp <- featuretable$data
  rm(featuretable)
  df.tmp <- as.data.frame(mat.tmp)
  rm(mat.tmp)
  df.tmp$OTUid <- rownames(df.tmp)
  rownames(df.tmp) <- NULL
  tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
}

#runfrom local: host.tmp <- tableimport.function("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.batASVs.table.qza")
host.tmp <- tableimport.function("host.qza")
colnames(host.tmp) <- c("ASVid", "SampleID", "Reads")

## import blast output (from `host_blastID.R` script):
## runfromlocal: tmp2 <- read_csv(file='~/Repos/mysosoup/data/host/Mangan.blastOut.csv.gz')
tmp2 <- read_csv(file = 'https://github.com/devonorourke/mysosoup/raw/master/data/host/Mangan.blastOut.csv')
host_species <- c('Myotis sodalis', 'Myotis lucifugus', 'Nycticeius humeralis')
tmp2 <- tmp2 %>% filter(species %in% host_species) %>% distinct(query, species)
tmp2 <- merge(host.tmp, tmp2, by.x='ASVid', by.y='query')

## get summary of read abundances per taxa for each sample
## summary of taxa counts:
SampTaxaReadSumry <- tmp2 %>% 
  group_by(SampleID, species) %>% 
  summarise(sumBatReads=sum(Reads)) %>%
  spread(., key = species, value = sumBatReads)

tmp1 <- tmp2 %>% 
  group_by(SampleID) %>% 
  summarise(sumBatReads=sum(Reads))

rm(tmp2, host.tmp)

## long workaround to flag instances with multiple species listed:
mylu <- is.na(SampTaxaReadSumry$`Myotis lucifugus`)
myso <- is.na(SampTaxaReadSumry$`Myotis sodalis`)
nyhu <- is.na(SampTaxaReadSumry$`Nycticeius humeralis`)

tmpout <- data.frame(SampTaxaReadSumry$SampleID, mylu, myso, nyhu)
rm(mylu, myso, nyhu)
colnames(tmpout)[1] <- "SampleID"
row.names(tmpout) <- tmpout$SampleID
tmpout$SampleID <- NULL
tmpout2 <- ifelse(tmpout == FALSE, 1, 0)
rm(tmpout)
tmpout2 <- as.data.frame(tmpout2)
tmpout2$rowsum <- rowSums(tmpout2[1:3])
tmpout2$SampleID <- row.names(tmpout2)
tmpout2 <- tmpout2 %>% select(rowsum, SampleID)
row.names(tmpout2 ) <- NULL
SampTaxaReadSumry <- merge(SampTaxaReadSumry, tmpout2)
rm(tmpout2)
SampTaxaReadSumry$multihit <- ifelse(SampTaxaReadSumry$rowsum == 1, FALSE, TRUE)
SampTaxaReadSumry <- SampTaxaReadSumry %>% select(-rowsum)
SampTaxaReadSumry <- merge(SampTaxaReadSumry, tmp1)
rm(tmp1)

## write to disk:
write.csv(SampTaxaReadSumry, file="~/Repos/mysosoup/data/taxonomy/host_abundances-by-taxaHits.csv", row.names = FALSE, quote = FALSE)

## add to summary file, include instances where no data exist:
# notrun: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan.nonbatASVs.table.qza", "nonbat.qza")
#runfrom local: nonbat.tmp <- tableimport.function("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.nonbatASVs.table.qza")
nonbat.tmp <- tableimport.function("arth.qza")
colnames(nonbat.tmp) <- c("ASVid", "SampleID", "Reads")
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv.gz", col_names = TRUE)
colnames(meta)[1] <- "SampleID"
df <- merge(nonbat.tmp, meta)
rm(nonbat.tmp, meta)
sumry0 <- as.data.frame(df %>% group_by(SampleID, SampleType, BatchType, Roost, CollectionMonth, Site) %>% 
                          summarise(sumNonBatReads=sum(Reads), nASVs=n_distinct(ASVid)) %>%
                          arrange(sumNonBatReads))
hostseq_data <- merge(sumry0, SampTaxaReadSumry, all=TRUE)

## write to disk:
write.csv(hostseq_data, file="~/Repos/mysosoup/data/taxonomy/sample_abundancesummaries_wBatHostData.csv", row.names = FALSE, quote = FALSE)

## summary of use
hostseq_data %>% 
  group_by(BatchType, multihit) %>% 
  summarise(counts=n()) %>% 
  spread(., key=multihit, value=counts)

bymeta <- hostseq_data %>% 
  group_by(BatchType, multihit, Roost, CollectionMonth, Site) %>% 
  summarise(counts=n()) %>%
  spread(., key=multihit, value=counts) %>%
  rename(BAT_multiHit_isFalse = `FALSE`, bat_multiHit_isTrue = `TRUE`, no_bat_data = `<NA>`) %>%
  arrange(BatchType, Site, Roost, CollectionMonth)
