library(tidyverse)
library(reshape2)
library(qiime2R)
library(scales)

# notrun: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan_noBats_ASVtable.qza", "arth.qza")

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


#runfrom local: dada2.tmp <- tableimport.function("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan_noBats_ASVtable.qza")
dada2.tmp <- tableimport.function("arth.qza")
colnames(dada2.tmp) <- c("ASVid", "SampleID", "Reads")

## add metadata  
#runfromlocal: meta <- read_csv(file = "~/Repos/mysosoup/data/manan_metadata.csv", col_names = TRUE)
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/manan_metadata.csv", col_names = TRUE)
colnames(meta)[1] <- "SampleID"

df <- merge(dada2.tmp, meta)
rm(dada2.tmp, meta)


## add taxonomy information
#runfromlocal: taxa <- read_delim(file = "~/Repos/mysosoup/data/taxonomy/mangan_tax_p97c94.tsv",delim = "\t", col_names = TRUE)
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv",
                   delim = "\t", col_names = TRUE)

taxa <- taxa %>% separate(., col = Taxon, sep=';',
                  into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>%
  select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
colnames(taxa)[1] <- "ASVid"

df <- merge(df, taxa)
rm(taxa)



### Filtering strategies
## distribution of number of ASVs per sample:
sumry0 <- df %>% group_by(SampleID, SampleType, BatchType) %>% summarise(sumReads=sum(Reads), nASVs=n_distinct(ASVid))

#1) drop singleton ASVs
## can create a filter that removes ASVs from dataset when they are present in only a single sample:
## number of ASVs per sample
sumry1 <- df %>% group_by(ASVid) %>% summarise(ASV_readCounts=sum(Reads), ASV_detections=n())
singletonASVs <- sumry1 %>% filter(ASV_detections > 1) %>% select(ASVid)
df_noSingleASV <- df %>% filter(ASVid %in% singletonASVs$ASVid)


#1b) retain only ASVs with at least Order-level taxonomy information
## also dropping the two ASVs assigned to a mule deer ... weird.
df_taxfiltd <- df %>% 
  filter(kingdom_name != "Unassigned") %>%
  filter(!is.na(phylum_name)) %>%
  filter(!is.na(class_name)) %>%
  filter(!is.na(order_name)) %>%
  filter(phylum_name != "Chordata")

dropASVs <- df %>% 
  filter(kingdom_name == "Unassigned" | phylum_name == "Chordata" | is.na(order_name)) %>%
  select(ASVid) %>% distinct(.)   ## 2213 ASVs 


#1c) also drop any read with less than N sequences (arbitrary)

df_taxfiltd_min10 <- df_taxfiltd %>% filter(Reads >= 10)
df_taxfiltd_min50 <- df_taxfiltd %>% filter(Reads >= 50)
df_taxfiltd_min100 <- df_taxfiltd %>% filter(Reads >= 100)

df_taxfiltd_min10 %>% select(SampleID) %>% n_distinct(.)  ## 296 samples
df_taxfiltd_min50 %>% select(SampleID) %>% n_distinct(.)  ## 296 samples
df_taxfiltd_min100 %>% select(SampleID) %>% n_distinct(.) ## 296 samples ... so a min-read filtering doesn't drop any particular sample

#2) same ideas as #1, but we group by NTC or not too:
sumry2 <- df %>% group_by(ASVid, SampleType) %>% summarise(ASV_readCounts=sum(Reads), ASV_detections=n())
#2b) same for tax-filtd df:
sumry2b <- df_taxfiltd %>% group_by(ASVid, SampleType) %>% summarise(ASV_readCounts=sum(Reads), ASV_detections=n())
#2c) same for abund-filt df:
sumry2ci <- df_taxfiltd_min10 %>% group_by(ASVid, SampleType) %>% summarise(ASV_readCounts=sum(Reads), ASV_detections=n())
sumry2cii <- df_taxfiltd_min50 %>% group_by(ASVid, SampleType) %>% summarise(ASV_readCounts=sum(Reads), ASV_detections=n())
sumry2ciii <- df_taxfiltd_min100 %>% group_by(ASVid, SampleType) %>% summarise(ASV_readCounts=sum(Reads), ASV_detections=n())

bigfilt_df <- df %>% filter(ASVid %in% sumry2ciii$ASVid)

#3) How many ASVs are private to control/sample subsets?
controlASVs <- sumry2 %>% filter(SampleType=="control") %>% distinct(ASVid)
sampleASVs <- sumry2 %>% filter(SampleType=="sample") %>% distinct(ASVid)
common <- intersect(controlASVs$ASVid, sampleASVs$ASVid)  ## 131 in common (so 20 of these ASVs are unique to controls)
#3b) Same for the tax-filtered df?
b.controlASVs <- sumry2b %>% filter(SampleType=="control") %>% distinct(ASVid)
b.sampleASVs <- sumry2b %>% filter(SampleType=="sample") %>% distinct(ASVid)
b.common <- intersect(b.controlASVs$ASVid, b.sampleASVs$ASVid)  ## 96 in common still! so most of what they're sharing are in fact highly likely to be arthropod samples, not weird contaminants


#4) Create data.frame of ONLY common ASVs...
df_common <- df %>% filter(ASVid %in% common)  
  ## makes up nearly 1/3 of all the data!
  ## can't just drop these ASVs completely...

df_common %>% group_by(SampleType) %>% summarise(ReadSums=sum(Reads), nSamples=n())

### Visualizing data
# 1) distribution of reads per ASV by Samples with ASV
ggplot(sumry1, aes(x=ASV_readCounts, y=ASV_detections)) + 
  geom_point()

# 2) distribution of reads per ASV by samples with ASV; highlighting NTC samples as separate from samples  
ggplot() + 
  geom_point(data=sumry2 %>% filter(SampleType=="sample"), 
             aes(x=ASV_readCounts, y=ASV_detections), alpha=0.4) +
  geom_point(data=sumry2 %>% filter(SampleType=="control"), 
             aes(x=ASV_readCounts, y=ASV_detections, color=SampleType), color='red', alpha=0.8) +
  scale_x_continuous(labels = comma)

# 3) same plot, but highlighting the ASVs shared between sample and controls (black are only in guano samples)
ggplot() + 
  geom_point(data=sumry2 %>% filter(!ASVid %in% b.common), 
             aes(x=ASV_readCounts, y=ASV_detections), alpha=0.4) +
  geom_point(data=sumry2 %>% filter(ASVid %in% b.common), 
             aes(x=ASV_readCounts, y=ASV_detections, color=SampleType), color='orange', alpha=0.8) +
  scale_x_continuous(labels=comma)

# 4) Histogram of read depth frequency among control and sample subsets 
ggplot(df_common, aes(x=Reads, fill=SampleType)) + 
  geom_histogram(position="dodge") +
  xlim(0,5000) ## change value as needed



##### plotting frequency of arthropod vs. chordate vs. unassigned; compare with whether or not they were flagged as reference chimeras


###



library(vegan)

## beta diversity function for RAREFIED data from FOX site for select MONTHS
rrarewithdrop <- function(x, sample) {
  rrarefy(x[rowSums(x) >= sample, , drop = FALSE], sample)
}

betadist.function <- function(betatest, binaryval) {
  tmp.meta <- bigfilt_df %>% distinct(SampleID)
  tmp.mat <- dcast(bigfilt_df, SampleID ~ ASVid, value.var = "Reads", fill=0)
  row.names(tmp.mat) <- tmp.mat$SampleID
  tmp.mat$SampleID <- NULL
  tmp.mat2 <- data.frame(rrarewithdrop(tmp.mat, 9313))  ## this drops out two guano samples
  rm(tmp.mat)
  namefilter <- row.names(tmp.mat2)
  i <- as.matrix(colSums(tmp.mat2, na.rm=T) != 0)
  tmp.mat3 <- tmp.mat2[, i]
  rm(tmp.mat2)
  tmp.meta <- tmp.meta %>% filter(SampleID %in% namefilter)  ## redo meta data to drop any samples that were discarded by rarefying
  tmp.betadiv <- vegdist(tmp.mat3, betatest, binary = binaryval)
}


dist <- betadist.function("bray", FALSE)
nmds_list <- metaMDS(dist, distance = "bray", trymax = 40, k = 4, autotransform = FALSE)
nmds_df <- as.data.frame(nmds_list$points)
rm(nmds_list)
nmds_df$SampleID <- row.names(nmds_df)

meta <- df %>% select(SampleID, CollectionMonth, Site, SampleType, BatchType, Roost) %>% distinct(.)

nmds_df <- merge(nmds_df, meta)


ggplot(nmds_df, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site)) + 
  geom_point() +
  facet_wrap(~SampleType)


nmds_df$CollectionMonth <- as.factor(nmds_df$CollectionMonth)
nmds_df$CollectionMonth <- as.factor(nmds_df$CollectionMonth)
