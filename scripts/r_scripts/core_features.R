## core features  -- identifying which ASVs are present in a high proportion of the group (or dataset)
library(tidyverse)
library(reshape2)
library(qiime2R)
library(scales)

################################################################################
## part 1 - loading read, meta, and taxa information
################################################################################

## import read data
# notrun: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan_noBats_ASVtable.qza", "nonbat.qza")

## function to convert .qza to matrix, then convert wide-format matrix to long-format data.frame object
tableimport.function <- function(table){
  featuretable <- read_qza(table)
  mat.tmp <- featuretable$data
  rm(featuretable)
  df.tmp <- as.data.frame(mat.tmp)
  rm(mat.tmp)
  df.tmp$OTUid <- rownames(df.tmp)
  rownames(df.tmp) <- NULL
  tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0) %>% 
  rename(ASVid = OTUid, SampleID = variable, Reads = value)
}

## run from local: tmp1 <- tableimport.function('~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.nonbatASVs.table.qza')
tmp1 <- tableimport.function('nonbat.qza')

## add metadata  
tmp2 <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv.gz", col_names = TRUE) %>% 
  select(SampleID, SampleType, CollectionMonth, Site) %>% 
  merge(tmp1, .) %>% 
  filter(SampleType == "sample") %>% 
  rename(Month = CollectionMonth)
tmp2$Site <- ifelse(tmp2$Site == "Egner", gsub("Egner", "EN", tmp2$Site), tmp2$Site)
tmp2$Site <- ifelse(tmp2$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", tmp2$Site), tmp2$Site)
tmp2$Month <- ifelse(tmp2$Month == "6", gsub("6", "June", tmp2$Month), tmp2$Month)
tmp2$Month <- ifelse(tmp2$Month == "7", gsub("7", "July", tmp2$Month), tmp2$Month)
tmp2$Month <- ifelse(tmp2$Month == "9", gsub("9", "September", tmp2$Month), tmp2$Month)
tmp2$SiteMonth <- paste(tmp2$Site, tmp2$Month, sep="-")

## add taxonomy information
tmp3 <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_vs.tsv.gz", delim = "\t", col_names = TRUE)
tmp3 <- tmp3 %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
tmp3 <- as.data.frame(apply(tmp3, 2, function(y) gsub(".__", "", y)))
tmp3 <- as.data.frame(apply(tmp3, 2, function(y) gsub("^$|^ $", NA, y)))
tmp3 <- as.data.frame(apply(tmp3, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
tmp3 <- as.data.frame(apply(tmp3, 2, function(y) gsub("Unassigned", NA, y)))
colnames(tmp3)[1] <- "ASVid"
data <- merge(tmp3, tmp2) %>% select(-SampleType)
rm(tmp1, tmp2, tmp3)


################################################################################
## part 2 - calculating shared features
################################################################################
tSamples <- n_distinct(data$SampleID)

coreASVsumry <- data.frame(data %>% group_by(ASVid, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>% 
  summarise(nSamples = n_distinct(SampleID), nReads = sum(Reads)) %>% 
  mutate(p5samp = (nSamples >= (round(0.1*tSamples)))) %>% 
  mutate(p10samp = (nSamples >= (round(0.1*tSamples)))) %>% 
  mutate(p20samp = (nSamples >= (round(0.2*tSamples)))) %>% 
  mutate(p30samp = (nSamples >= (round(0.3*tSamples)))) %>% 
  mutate(p40samp = (nSamples >= (round(0.4*tSamples)))) %>% 
  mutate(p50samp = (nSamples >= (round(0.5*tSamples)))) %>% 
  mutate(p60samp = (nSamples >= (round(0.6*tSamples)))) %>% 
  mutate(p70samp = (nSamples >= (round(0.7*tSamples)))) %>% 
  mutate(p80samp = (nSamples >= (round(0.8*tSamples)))))

## which of these are not what we expect?
p10samps <- coreASVsumry %>% filter(p10samp == TRUE)
p10oddballs <- p10samps %>% filter(is.na(phylum_name) | phylum_name != "Arthropoda")
p10oddASVs <- data.frame(p10oddballs) %>% select(ASVid) %>% mutate_if(is.factor, as.character) %>% pull()
  ## manually BLASTed these 9 ASVs and only one is a likely arthropod; 
  ## others either undetectable or have pid < 0.8, or are fungus/mold
  ## going to restrict our analyses to those sequences with matching Arthropod taxa and family-level info

filt_samps <- coreASVsumry %>% filter(phylum_name == "Arthropoda") %>% filter(!is.na(family_name))

## how many of each taxa Order are represented in the Nth percentile?
perc5tbl <- filt_samps %>% filter(p5samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p05')
perc10tbl <- filt_samps %>% filter(p10samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p10')
perc20tbl <- filt_samps %>% filter(p20samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p20')
perc30tbl <- filt_samps %>% filter(p30samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p30')
perc50tbl <- filt_samps %>% filter(p50samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p50')
perc_alltbl <- rbind(perc5tbl, perc10tbl, perc20tbl, perc30tbl, perc50tbl)

perc_alltbl <- perc_alltbl %>% spread(order_name, counts, fill=NA) %>% arrange(percentile)
write.table(perc_alltbl, 
            file="~/Repos/mysosoup/data/text_tables/core_features/percentile_corefeatures_table.txt",
            quote=FALSE, row.names = FALSE)

## which ASVs and taxa are among these top 5th percentile features?
perc5names <- filt_samps %>% 
  filter(p5samp == TRUE)
write.table(perc5names, 
            file="~/Repos/mysosoup/data/text_tables/core_features/summary_corefeatures_table.txt",
            quote=FALSE, row.names=FALSE)
