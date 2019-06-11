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
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "6", gsub("6", "June", tinymeta$CollectionMonth), tinymeta$CollectionMonth)
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "7", gsub("7", "July", tinymeta$CollectionMonth), tinymeta$CollectionMonth)
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "9", gsub("9", "September", tinymeta$CollectionMonth), tinymeta$CollectionMonth)

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

## Beta diversity
## calculate distances and generate one plot per PID dataset
## Facet the 3 axes per NMDS plot, per Distance method
NMDSplotfunction <- function(physeqdata, DistTest, BinaryTest, DistMethod) {
  rphy_mat = as(otu_table(physeqdata), "matrix")
  dist_out <- vegdist(x = rphy_mat, method = DistTest, binary = BinaryTest)
  nmds_out <- metaMDS(comm = dist_out, distance = DistTest, k = 3, trymax = 30)
  data.frame(nmds_out$points) %>% mutate(SampleID = row.names(.)) %>% mutate(DistMethod=DistMethod)
}

#### ---- RESET THESE VALUES FOR EACH PID ---- ####
## ex. {nmds_df.dice <- NMDSplotfunction(rphy_p99, "bray", TRUE, "dice")}
nmds_df.dice <- NMDSplotfunction(rphy_p100, "bray", TRUE, "dice")
nmds_df.bray <- NMDSplotfunction(rphy_p100, "bray", FALSE, "bray")
nmds_df.mori <- NMDSplotfunction(rphy_p100, "morisita", FALSE, "morisita")
#### ---- RESET ABOVE VALUES FOR EACH PID ---- ####

nmds_df.all <- rbind(nmds_df.dice, nmds_df.bray, nmds_df.mori)
nmds_df.all <- merge(nmds_df.all, tinymeta)

nmds_df.all$DistMethod <- factor(nmds_df.all$DistMethod, levels = c("dice", "bray", "morisita"))
nmds_df.all$CollectionMonth <- factor(nmds_df.all$CollectionMonth, levels = c("June", "July", "September", "control"))

v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)
p.nmds12 <- ggplot(nmds_df.all, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site)) + 
  geom_point() +
  facet_grid(~ DistMethod) +
  scale_color_manual(values = c(v3pal, "gray60")) +
  scale_x_continuous(limits = c(-.55, .65)) +
  scale_y_continuous(limits = c(-.55, .65)) +
  theme_devon() +
  theme(legend.position = "none")
p.nmds13 <- ggplot(nmds_df.all, aes(x=MDS1, y=MDS3, color=CollectionMonth, shape=Site)) + 
  geom_point() +
  facet_grid(~ DistMethod) +
  scale_color_manual(values = c(v3pal, "gray60")) +
  scale_x_continuous(limits = c(-.55, .65)) +
  scale_y_continuous(limits = c(-.55, .65)) +
  theme_devon() +
  theme(legend.position = "none", strip.background.x = element_blank(), strip.text.x = element_blank())
p.nmds23 <- ggplot(nmds_df.all, aes(x=MDS2, y=MDS3, color=CollectionMonth, shape=Site)) + 
  geom_point() +
  facet_grid(~ DistMethod) +
  scale_color_manual(values = c(v3pal, "gray60")) +
  scale_x_continuous(limits = c(-.55, .65)) +
  scale_y_continuous(limits = c(-.55, .65)) +
  theme_devon() +
  theme(legend.position = "none", strip.background.x = element_blank(), strip.text.x = element_blank())

p.fake <- ggplot(nmds_df.all, aes(x=MDS1, y=MDS3, color=CollectionMonth, shape=Site)) + 
  geom_point() +
  facet_grid(~ DistMethod) +
  scale_color_manual(values = c(v3pal, "gray60")) +
  theme_devon()

legend_b <- get_legend(p.fake + theme(legend.position="bottom"))

plot_grid(p.nmds12, p.nmds13, p.nmds23, legend_b, nrow = 4, rel_heights = c(1.1, 1, 1, 0.2))
