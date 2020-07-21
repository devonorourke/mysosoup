library(tidyverse)
library(reshape2)
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



# notrun: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan_noBats_ASVtable.qza", "nonbat.qza")

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


## runfrom local: nonbat.tmp <- tableimport.function("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.nonbatASVs.table.qza")
nonbat.tmp <- tableimport.function("nonbat.qza")
colnames(nonbat.tmp) <- c("ASVid", "SampleID", "Reads")

## add metadata
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE)
df <- merge(nonbat.tmp, meta)
rm(nonbat.tmp, meta)

## add taxonomy information
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"
df <- merge(df, taxa)

## calculate the total number of reads per sample, and the number of unique ASVs observed per sample
sumry0 <- df %>% 
  group_by(SampleID, SampleType, BatchType, Source, ContamArea) %>% 
  summarise(sumReads=sum(Reads), nASVs=n_distinct(ASVid))
write.csv(sumry0, file = "~/Repos/mysosoup/data/text_tables/contam_evals/freqObs_and_seqDepth_perSample.csv",
          quote = FALSE, row.names = FALSE)
## 5 of 7 negative control (NTC) samples have very low number of ASVs but fairly average numbers of reads per sample... most guano have >30 ASVs
## only 2 of 7 true NTCs has more than 30, just like the "blankS39" sample which we KNOW had guano fall into the well

## related question: how many times do we see a particular ASV among any of the NTC samples?
## are those NTCs with few ASVs at least all showing the same ASVs?
df %>% filter(SampleType == "control") %>% summarise(nSamples = n_distinct(SampleID))   ## there are 8 NTC samples (7 extraction blanks, one we know is contaminated ("blankS39")...

NTConly_df <- df %>% 
  filter(SampleType == "control") %>%
  group_by(ASVid) %>% 
  summarise(counts=n())
## nope; we find that only 4/87 ASVs in four of eight samples
## 5/87 ASVs in three samples; 15/87 in two samples
## 63/87 in only single sample ... so most are unique, indicating that there isn't pervasive (ex. reagent) contamination


## Coloring the samples that were withing a cell of an NTC
## any relationship between the number of ASVs per sample, the number of reads per sample, and whether or not
## ..the sample was derived from a Plate-extraction?

## plot this; save as 'contam_eval_by_PlateType_and_ContamZone'; export at 
ggplot(sumry0, aes(sumReads, nASVs, color=ContamArea, shape=SampleType)) + 
  geom_point(data=sumry0 %>% filter(SampleType == "control"), aes(sumReads, nASVs, color=ContamArea, shape=SampleType), size=4) + 
  geom_point(data=sumry0 %>% filter(SampleType != "control"), aes(sumReads, nASVs, color=ContamArea, shape=SampleType)) + 
  facet_wrap(~ Source) +
  scale_x_continuous(labels = comma) +
  theme_devon()
## these results indicate that the kind of extraction (plate vs. isolate) isn't really diving differences in #ASVs or #reads per sample
## location of true sample near a contaminated well doesn't influence number of reads or number of ASVs - they're indendent



## drop singleton ASVs
## can create a filter that removes ASVs from dataset when they are present in only a single sample:
## number of ASVs per sample
asv_ReadCounts <- df %>% group_by(ASVid) %>%  summarise(ASV_readCounts=sum(Reads), ASV_detections=n())
singletonASVs <- asv_ReadCounts %>% filter(ASV_detections == 1) %>% select(ASVid)
df_noSingleASV <- df %>% filter(!ASVid %in% singletonASVs$ASVid)

## make same plot as above: do we see similar trends with plate vs. isolates and controls vs true samples when singletons are removed?
sumry1 <- df_noSingleASV %>% 
  group_by(SampleID, SampleType, BatchType, Source, ContamArea) %>% 
  summarise(sumReads=sum(Reads), nASVs=n_distinct(ASVid))

ggplot(sumry1, aes(sumReads, nASVs, color=ContamArea, shape=SampleType)) + 
  geom_point(data=sumry1 %>% filter(SampleType == "control"), aes(sumReads, nASVs, color=ContamArea, shape=SampleType), size=4) + 
  geom_point(data=sumry1 %>% filter(SampleType != "control"), aes(sumReads, nASVs, color=ContamArea, shape=SampleType)) + 
  facet_wrap(~ Source) +
  scale_x_continuous(labels = comma) +
  theme_devon()
## results are similar to previous plot... filtering out singletons doesn't change overall patterns

## drop any ASV that doesn't have at least Family-rank information (assumption here is poor classification MAY be indicative of sequence error)
df_taxfilt <- df %>% 
  filter(!is.na(kingdom_name)) %>% 
  filter(!is.na(phylum_name)) %>% 
  filter(!is.na(class_name)) %>% 
  filter(!is.na(order_name)) %>% 
  filter(!is.na(family_name)) %>% 
  filter(phylum_name == "Arthropoda")
sumry_taxfilt <- df_taxfilt %>% group_by(ASVid) %>% summarise(sumReads=sum(Reads), nSamples=n_distinct(SampleID))
commonASVs_taxfilt <- intersect(df_taxfilt %>% filter(SampleType == "control") %>% select(ASVid),
                                df_taxfilt %>% filter(SampleType == "sample") %>% select(ASVid)) %>% pull()

unique_controlASVs_taxfilt <- setdiff(df_taxfilt %>% filter(SampleType == "control") %>% select(ASVid),
                                      df_taxfilt %>% filter(SampleType == "sample") %>% select(ASVid)) %>% pull()
unique_sampleASVs_taxfilt <- setdiff(df_taxfilt %>% filter(SampleType == "sample") %>% select(ASVid),
                                     df_taxfilt %>% filter(SampleType == "control") %>% select(ASVid)) %>% pull()

## plot; save as 'contam_eval_ASV_and_SeqCounts_uniqnessColored_per_Guano_or_Control'
ggplot(sumry_taxfilt, aes(x=sumReads, y=nSamples)) +
  geom_point(data = sumry_taxfilt %>% filter(ASVid %in% commonASVs_taxfilt), color="black") +
  geom_point(data = sumry_taxfilt %>% filter(ASVid %in% unique_sampleASVs_taxfilt), color="blue") +
  geom_point(data = sumry_taxfilt %>% filter(ASVid %in% unique_controlASVs_taxfilt), color="red", size=2) +
  theme_devon()


## what's up with those common ASVs in both true samples and contaminants?
common_ASV_df <- df %>% filter(ASVid %in% commonASVs_taxfilt)
sumry_common_contamASVs <- common_ASV_df %>% group_by(ASVid) %>% summarise(nReads=sum(Reads), nSamples=n_distinct(SampleID))
## how many total sequences are dropped if we remove those 69 ASVs?
df %>% summarise(nReads=sum(Reads))
common_ASV_df %>% summarise(nReads=sum(Reads))
4253793/9411633 ## these 68 'contaminant' ASVs represent nearly half (45.2%) of all ASVs !


##### filtering options
## filter by 1) removing ASVs that are identified in more NTCs than true samples (by occurrence, not percentage) .. and/OR ..
##           2) remove ASVs that contain more total reads for that ASV among all NTCs than among all true samples

contam_df <- df %>% filter(SampleType == "control")

noncontam_sumry <- df %>% 
  filter(SampleType == "sample") %>% 
  group_by(ASVid) %>% 
  summarise(NC_nReads=sum(Reads), NC_nSamples=n_distinct(SampleID))

contam_sumry <- contam_df %>% group_by(ASVid) %>% 
  summarise(C_nReads=sum(Reads), C_nSamples=n_distinct(SampleID))

light_filt_df <- merge(contam_sumry, noncontam_sumry)
light_filt_df$more_ASVs <- light_filt_df$C_nSamples > light_filt_df$NC_nSamples  ## ZERO examples
light_filt_df$more_Reads <- light_filt_df$C_nReads > light_filt_df$NC_nReads  ## Just 9 examples here

## these filtering strategies suggest that despite some ASVs occurring in many samples, the ASVs that are shared between ..
## .. contaminant and true samples rarely produce what you'd expect as indicative of widespread contamination: all NTCs having the most widely observed ASVs
## Instead, we find that there are NEVER more than 4 NTCs with these common ASVs (just 4 of 79 ASVs); five examples where there are 3 NTCs with a common ASV; ..
## .. and the remaining 70 ASVs are only ever detected in 2 or 1 NTC samples.

## we can further refind this by highlight just the ASVs which would be used in downstream analyses for diet work: ..
## .. those ASVs that contain Arthropod phylum-classification, and at least Family-rank information:

taxfilt_contamASVs <- intersect(light_filt_df$ASVid, df_taxfilt$ASVid)
contam_taxfilt_df <- light_filt_df %>% filter(ASVid %in% taxfilt_contamASVs)

## this further reduces our analyses to just 68 ASVs

## We'll create two lists of ASVs to filter our original dataset:
## 1) ASVs that meet our taxonomy completeness criteria (Arthropod-associated phylum_name with at least Family info) and include all ASVs whether or not they are present in NTC samples
## 2) As with #1 above, but removing the 68 ASVs present in NTC samples 

taxfilt_ASVs <- df_taxfilt %>% distinct(ASVid)
colnames(taxfilt_ASVs)[1] <- "#OTUID"
write.table(taxfilt_ASVs, "~/github/mysosoup/data/taxonomy/taxfiltd_ASVs_NTCincluded.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)

taxfilt_2_ASVs <- df_taxfilt %>% filter(!ASVid %in% contam_df$ASVid) %>% distinct(ASVid)
colnames(taxfilt_2_ASVs)[1] <- "#OTUID"
write.table(taxfilt_2_ASVs, "~/github/mysosoup/data/taxonomy/taxfiltd_ASVs_NTCdrops.txt", row.names = FALSE, quote = FALSE, col.names = TRUE)

## These two text files are used to filter by retaining only the ASVs listed in those files. 
## The resulting `.qza` artifacts are then rarefied and then used for diversity analyses.

## double check to ensure there are NTCs (or not) in these two filtered ASV sets:
df %>% filter(ASVid %in% taxfilt_ASVs$ASVid)
df %>% filter(ASVid %in% taxfilt_2_ASVs$ASVid)
