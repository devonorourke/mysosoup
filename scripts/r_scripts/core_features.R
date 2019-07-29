## core features  -- identifying which ASVs are present in a high proportion of the group (or dataset)
library(tidyverse)
library(reshape2)
library(qiime2R)
library(scales)
library(formattable)
library(htmltools)
library(webshot)

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

## need ASValias too...
tmpqzaimport.function <- function(qzapath){
  featuretable <- read_qza(qzapath)
  mat.tmp <- featuretable$data
  rm(featuretable)
  df.tmp <- as.data.frame(mat.tmp)
  rm(mat.tmp)
  df.tmp$OTUid <- rownames(df.tmp)
  rownames(df.tmp) <- NULL
  tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
  colnames(tmp) <- c("ASVid", "SampleID", "Reads")
  tmp
}
## download file: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/seqs/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza", "nonfyd.qza")
# save path for each file.. for example: nonpath="~/Downloads/nonfyd.qza"
## or run from local: nonpath="/Users/do/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza"
nonfyd_reads <- tmpqzaimport.function(nonpath)
tmp4 <- nonfyd_reads %>% group_by(ASVid) %>%  summarise(nReads=sum(Reads)) %>% arrange(-nReads) %>% mutate(ASValias=paste0("ASV-", row.names(.))) %>% select(-nReads)
rm(nonfyd_reads)
data <- merge(data, tmp4, all.x = TRUE) %>% 
  drop_na(ASValias) ## this drops 3 rows where a couple of unrarefied samples with just one ASV per sample were dropped... 
rm(tmp4)

################################################################################
## part 2 - calculating shared features
################################################################################
tSamples <- n_distinct(data$SampleID)

coreASVsumry <- data.frame(data %>% group_by(ASValias, phylum_name, class_name, order_name, family_name, genus_name, species_name) %>% 
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
perc5names <- coreASVsumry %>% 
  filter(p5samp == TRUE)
write.table(perc5names, 
            file="~/Repos/mysosoup/data/text_tables/core_features/summary_corefeatures_table.txt",
            quote=FALSE, row.names=FALSE)

## which Orders have the most of the ASVs we see?
perc5names %>% 
  group_by(order_name) %>% 
  summarise(counts=n_distinct(ASValias)) %>% 
  arrange(-counts)
  ## 44 of the 63 are Dipteran!
  
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
  select(ASValias, order_name, family_name, genus_name, species_name, nSamples, nReads) %>% 
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
  select(-ASValias) %>% 
  arrange(-nSamples) %>% 
  mutate(Species = gsub(" ", "_", .$Species))
write.table(perc5names_reduced_noASV, 
            file="~/Repos/mysosoup/data/text_tables/core_features/corefeatures_simple_noASV.txt",
            quote=FALSE, row.names=FALSE)

## try setting up a table with top 10% targets (same ASVs as top 5%), but add in 30, 50, 70%, and highlight with red/green FALSE/TRUE in columns
prettyTable <- coreASVsumry %>% 
  filter(p20samp == TRUE) %>% 
  select(ASValias, order_name, family_name, genus_name, species_name, nSamples, nReads, p30samp, p50samp, p70samp) %>% 
  rename(Order=order_name, Family=family_name, Genus=genus_name, Species=species_name,
         #Top30p=p30samp, Top50p=p50samp, Top70p=p70samp,
         `Top 30 %`=p30samp, `Top 50 %`=p50samp, `Top 70 %`=p70samp,
         Samples=nSamples, SeqCounts=nReads) %>% 
  arrange(-`Top 70 %`, -`Top 50 %`, -`Top 30 %`, Order, Family, Genus, Species, Samples)
prettyTable$Order <- as.character(prettyTable$Order)
prettyTable$Family <- as.character(prettyTable$Family)
prettyTable$Genus <- as.character(prettyTable$Genus)
prettyTable$Species <- as.character(prettyTable$Species)
prettyTable[is.na(prettyTable)] <- ""

## function to export formattable object properly
export_formattable <- function(f, file, width = "95%", height = NULL,
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}


## export here; save as 'CoreFeautres'
setwd("~/Repos/mysosoup/figures/")

corefeat <- formattable(prettyTable, list(
  Samples = color_tile("white", "#ff796c"),
  SeqCounts = color_tile("white", "#bf9005"),
  Top30p = formatter("span",
                         style = x ~ style(color = ifelse(x, "green", "red")),
                         x ~ icontext(ifelse(x, TRUE, FALSE), ifelse(x, "Yes", "No"))),
  Top50p = formatter("span",
                     style = x ~ style(color = ifelse(x, "green", "red")),
                     x ~ icontext(ifelse(x, TRUE, FALSE), ifelse(x, "Yes", "No"))),
  Top70p = formatter("span",
                     style = x ~ style(color = ifelse(x, "green", "red")),
                     x ~ icontext(ifelse(x, TRUE, FALSE), ifelse(x, "Yes", "No")))
  
))


export_formattable(corefeat,"CoreFeautres.png")
