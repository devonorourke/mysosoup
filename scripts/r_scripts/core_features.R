## core features  -- identifying which ASVs are present in a high proportion of the group (or dataset)
library(tidyverse)
library(reshape2)
library(qiime2R)
library(scales)
library(formattable)

################################################################################
## part 1 - loading read, meta, and taxa information
################################################################################

## import read data
# notrun: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza", "nonbat.qza")

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

## run from local: tmp1 <- tableimport.function('~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza')
tmp1 <- tableimport.function('nonbat.qza')

## add metadata  
tmp2 <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE) %>% 
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
tmp3 <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
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

## how many of each taxa Order are represented in at least 5%, 10% 50% of samples?
perc5tbl <- p10samps %>% filter(p5samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p05')
perc10tbl <- p10samps %>% filter(p10samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p10')
perc20tbl <- p10samps %>% filter(p20samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p20')
perc30tbl <- p10samps %>% filter(p30samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p30')
perc50tbl <- p10samps %>% filter(p50samp == TRUE) %>% group_by(order_name) %>% summarise(counts=n()) %>% mutate(percentile='p50')
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

## which Orders have the most of the ASVs we see?
perc5names %>% 
  group_by(order_name) %>% 
  summarise(counts=n_distinct(ASVid)) %>% 
  arrange(-counts)

## which Orders have the most of the Families we see?
perc5names %>% 
  group_by(order_name) %>% 
  summarise(counts=n_distinct(family_name)) %>% 
  arrange(-counts)

## which Families are most frequent?
perc5names %>% 
  group_by(order_name, family_name) %>% 
  summarise(counts=n()) %>% 
  arrange(-counts)

## reduce the `perc5names` dataset to just taxa info and nSample/nRead data:
perc5names_reduced <- perc5names %>% 
  select(ASVid, order_name, family_name, genus_name, species_name, nSamples, nReads) %>% 
  rename(Order=order_name, Family=family_name, Genus=genus_name, Species=species_name) %>% 
  arrange(-nSamples)

perc5names_reduced <- perc5names_reduced %>% 
  mutate_if(is.factor, as.character)

perc5names_reduced[is.na(perc5names_reduced)] <- ""

## export raw text table:
write.table(perc5names_reduced, 
            file="~/Repos/mysosoup/data/text_tables/core_features/corefeatures_simple.txt",
            quote=FALSE, row.names=FALSE)

## or save as fancy colored table/plot;
## could save as web page, but don't have current infrastructure to embed...
formattable(perc5names_reduced, list(
  nSamples = color_tile("white", "#ff796c"),
  nReads = color_tile("white", "#bf9005")
))


## export raw text table without any ASV info (to fit on single page):
perc5names_reduced_noASV <- perc5names_reduced %>% 
  select(-ASVid) %>% 
  arrange(-nSamples) %>% 
  mutate(Species = gsub(" ", "_", .$Species))
write.table(perc5names_reduced_noASV, 
            file="~/Repos/mysosoup/data/text_tables/core_features/corefeatures_simple_noASV.txt",
            quote=FALSE, row.names=FALSE)