library(qiime2R)
library(tidyverse)
library(reshape2)
library(formattable)


## code extended from `machine_learn_analyses.R` script
## import read datasets for rarefied data:
qzaimport.function <- function(qzapath){
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

## download file: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/seqs/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza", "rarfyd.qza")
# save path for each file.. for example: 
#    rfypath="~/Downloads/rarfyd.qza"
# alternative: run from local rfypath="/Users/do/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza"
rarfyd_reads <- qzaimport.function(rfypath)

## add taxonomy information
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"
row.names(taxa) <- taxa$ASVid

## import metadata
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE) %>% 
  filter(SampleType == "sample") %>% select(SampleID, BatchType, Roost, CollectionMonth, Site)
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "6", gsub("6", "June", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "7", gsub("7", "July", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "9", gsub("9", "September", meta$CollectionMonth), meta$CollectionMonth)
meta$SiteMonth <- paste(meta$Site, meta$CollectionMonth, sep="-")

## merge the data
rfy_plotdat <- merge(rarfyd_reads, taxa) %>% merge(., meta)
rm(meta, taxa, rarfyd_reads, rfypath, qzaimport.function)


## summary describes the number of detections and read abundances per ASV
## however, these are not aggregated per Sample ...
ordersum1 <- rfy_plotdat %>% 
  group_by(order_name) %>% 
  summarise(OrderReads=sum(Reads), OrderSamples=n()) %>% 
  mutate(pOrderSamples=round(OrderSamples/sum(OrderSamples), 3), pOrderReads=round(OrderReads/sum(OrderReads), 3)) %>% 
  arrange(-OrderSamples)

## This table depicts the perSample basis.
tmp1 <- rfy_plotdat %>% 
  group_by(SampleID, order_name) %>% 
  summarise(nOrder=n()) %>%
  spread(key=order_name, value=nOrder)
row.names(tmp1) <- tmp1$SampleID
tmp1$SampleID <- NULL
tmp1[!is.na(tmp1)] <- 1
tmp1$SampleID <- row.names(tmp1)
row.names(tmp1) <- NULL

tmp2 <- tmp1 %>% 
  melt(id.vars="SampleID") %>% 
  group_by(variable) %>%
  filter(!is.na(value)) %>%
  summarise(Samples=sum(value)) %>% 
  mutate(fracSamples=round(Samples/nrow(tmp1), 3)) %>% 
  arrange(-fracSamples)

## write to disk
write.csv(tmp2, 
          file="~/Repos/mysosoup/data/text_tables/taxa_summaries/perOrder_Sample_andRead_sumry.csv",
          quote=FALSE, row.names = FALSE)

selectOrders <- tmp2 %>% filter(fracSamples > 0.05) %>% select(variable) %>% mutate_if(is.factor, as.character) %>% pull()
tmp2select <- tmp2 %>% filter(fracSamples > 0.05) %>% rename(order_name = variable)

## how many distinct families
rfy_plotdat %>% 
  filter(order_name %in% selectOrders) %>% 
  group_by(family_name) %>% 
  summarise(nFamilies=n()) %>% nrow()
  ## 190

tmp3 <- rfy_plotdat %>% 
  filter(order_name %in% selectOrders) %>% 
  group_by(order_name) %>% 
  summarise(Families=n_distinct(family_name), Genera=n_distinct(genus_name), 
            Species=n_distinct(species_name), ASVs=n_distinct(ASVid)) %>% 
  arrange(-Families)

tmp4 <- merge(tmp2select, tmp3) %>% select(-fracSamples) %>% arrange(-Samples) %>% rename(`Arthropod Order`=order_name)

## write to disk
write.csv(tmp4, 
          file="~/Repos/mysosoup/data/text_tables/taxa_summaries/perSelectOrder_Fam-Gen-Spec_sumry.csv",
          quote=FALSE, row.names = FALSE)

## make pretty with formattable package and export as plot instead of data table;
## plot; save as 'perSelectOrder_Fam-Gen-Spec-ASV_sumry'; export at 700x550
formattable(tmp4, list(
            Families = color_tile("white", "orange"),
            Genera = color_tile("white", "violet"),
            Species = color_tile("white", "dodgerblue"),
            ASVs = color_tile("white", "seagreen3")))
