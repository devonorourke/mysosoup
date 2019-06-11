library(tidyverse)
library(vegan)
library(phyloseq)
library(cowplot)
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

## import metadata 
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/mangan_metadata.csv.gz", col_names = TRUE)
tinymeta <- meta %>% select(SampleID, Roost, CollectionMonth, Site)
tinymeta$Site <- ifelse(tinymeta$Site == "Egner", gsub("Egner", "EN", tinymeta$Site), tinymeta$Site)
tinymeta$Site <- ifelse(tinymeta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", tinymeta$Site), tinymeta$Site)
tinymeta$CollectionMonth[is.na(tinymeta$CollectionMonth)] <- "control"
tinymeta$Labeler <- paste(tinymeta$Site, tinymeta$Roost, sep="-")
tinymeta$Labeler <- ifelse(tinymeta$Labeler == "control-control", gsub("control-control", "control", tinymeta$Labeler), tinymeta$Labeler)
tinymeta$CollectionMonth <- as.factor(tinymeta$CollectionMonth)


## import taxonomy info
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv.gz", 
                   col_names = TRUE, delim = "\t")
taxa <- taxa %>% separate(., col = Taxon, sep=';',
                          into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>%
  select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
colnames(taxa)[1] <- "ASVid"


## import sequence data:
## import 3 datasets: data clustered with 100% identity, 99%, and 98%
## function filters out any ASV that isn't an Arthropod with at least Family-rank information

tableimport.function <- function(table){
  featuretable <- read_qza(table)
  mat.tmp <- featuretable$data
  rm(featuretable)
  df.tmp <- as.data.frame(mat.tmp)
  rm(mat.tmp)
  df.tmp$OTUid <- rownames(df.tmp)
  rownames(df.tmp) <- NULL
  tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
  colnames(tmp) <- c("ASVid", "SampleID", "Reads")
  tmp %>% filter(ASVid %in% (taxa %>% filter(phylum_name == "Arthropoda") %>% filter(!is.na(family_name)) %>% select(ASVid) %>% pull()))
}

dat_p100 <- tableimport.function("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan_noBats_ASVtable.qza")
dat_p99 <- tableimport.function("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan_noBats_OTUtable_pid99.qza")
dat_p98 <- tableimport.function("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan_noBats_OTUtable_pid98.qza")

## create physeq metadata and taxonomy data for analyses:
row.names(meta) <- meta$SampleID
phy_meta <- sample_data(meta)

row.names(taxa) <- taxa$ASVid
tax_mat <- as.matrix(taxa)
phy_tax <- tax_table(tax_mat)
rm(tax_mat)

## Create unique physeq OTUtable objects for each dataset..
## and combine with Meta/Taxa data to create single physeq object.
physeqImport <- function(data){
  data_mat <- dcast(data, SampleID ~ ASVid, value.var = 'Reads', fill = 0)  
  row.names(data_mat) <- data_mat$SampleID
  NULL -> data_mat$SampleID
  OTU <- otu_table(data_mat, taxa_are_rows = FALSE) ## import as physeq object 
  phyloseq(OTU, phy_tax, phy_meta)
}

phy_p100 <- physeqImport(dat_p100)
phy_p99 <- physeqImport(dat_p99)
phy_p98 <- physeqImport(dat_p98)

rm(phy_tax, phy_meta)

## rarefy:
## to estimate sampling depth, evaluate number of raw reads per sample:
p100_Sumry <- dat_p100 %>% group_by(SampleID) %>% summarise(nASVs=n_distinct(ASVid), nReads=sum(Reads))
p99_Sumry <- dat_p99 %>% group_by(SampleID) %>% summarise(nASVs=n_distinct(ASVid), nReads=sum(Reads))
p98_Sumry <- dat_p98 %>% group_by(SampleID) %>% summarise(nASVs=n_distinct(ASVid), nReads=sum(Reads))

## visualize distribution
p <- ggplot(p100_Sumry, aes(nReads)) + geom_histogram(fill="gray20", color="red", bins=20) + scale_x_continuous(labels = comma) + xlim(0,20000)
q <- ggplot(p99_Sumry, aes(nReads)) + geom_histogram(fill="gray50", color="red", bins=20) + scale_x_continuous(labels = comma) + xlim(0,20000)
r <- ggplot(p98_Sumry, aes(nReads)) + geom_histogram(fill="gray80", color="red", bins=20) + scale_x_continuous(labels = comma) + xlim(0,20000)
plot_grid(p, q, r, nrow = 3)

## arbitrarily picking 9000 because it appears to capture most data without going into left tail of distribution
rphy_p100 = rarefy_even_depth(phy_p100, sample.size=9000, replace = FALSE, rngseed = 1000)
rphy_p99 = rarefy_even_depth(phy_p99, sample.size=9000, replace = FALSE, rngseed = 1000)
rphy_p98 = rarefy_even_depth(phy_p98, sample.size=9000, replace = FALSE, rngseed = 1000)

nsamples(rphy_p100)  ## 276 samples remain (from 296)
nsamples(rphy_p99)  ## 276 samples remain (from 295)
nsamples(rphy_p98)  ## 276 samples remain (from 295)
  ## So, all 3 datasets have same number of samples

ntaxa(rphy_p100) ## 2871 taxa remain (from 2988)
ntaxa(rphy_p99) ## 1316 taxa remain (from 2988)
ntaxa(rphy_p98) ## 1031 taxa remain (from 2988)
  ## 3 datasets have very differet number of ASVs

## Alpha diversity
## Export rarefied Physeq OTU table back as matrix..
## then estimate Hill Numbers for q=0, q=1, and q=2
HillCalculator <- function(physeqdata, pid) {
  rphy_mat = as(otu_table(physeqdata), "matrix")  
  tmp.hill <- data.frame(renyi(rphy_mat, scales = c(0,1,2), hill=TRUE)) %>%
    mutate(SampleID = row.names(.))
  colnames(tmp.hill)[1:3] <- c("q=0", "q=1", "q=2")
  Alpha_df <- gather(tmp.hill, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2')) %>%
    mutate(clustPID = pid)
  merge(Alpha_df, tinymeta)
}

hill_p100 <- HillCalculator(rphy_p100, "pid-100")
hill_p99 <- HillCalculator(rphy_p99, "pid-99")
hill_p98 <- HillCalculator(rphy_p98, "pid-98")

hill_all <- rbind(hill_p100, hill_p99, hill_p98)
rm(hill_p100, hill_p99, hill_p98)

hill_all$CollectionMonth <- factor(hill_all$CollectionMonth, levels = c("June", "July", "September", "control"))
hill_all$clustPID <- factor(hill_all$clustPID, levels = c('pid-100', 'pid-99', 'pid-98'))

v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)

## plot alpha diversity
ggplot(hill_all, aes(x=Labeler, y=Hill_value, color=CollectionMonth)) + 
  geom_jitter(width = 0.1, alpha=0.8) + 
  scale_color_manual(values = c(v3pal, "gray40")) +
  facet_grid(clustPID ~ Hill_qType) +
  labs(x="", y="Estimated diversity") +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1),
        legend.position = "top")

