library(tidyverse)
library(phyloseq)
library(reshape2)
library(vegan)
library(ape)


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

## arbitrarily picking 9000 because it appears to capture most data without going into left tail of distribution
rphy_p100 = rarefy_even_depth(phy_p100, sample.size=9000, replace = FALSE, rngseed = 1000)
rphy_p99 = rarefy_even_depth(phy_p99, sample.size=9000, replace = FALSE, rngseed = 1000)
rphy_p98 = rarefy_even_depth(phy_p98, sample.size=9000, replace = FALSE, rngseed = 1000)

rm(phy_p100, phy_p99, phy_p98)

## set up 2 functions to calculate distances:
## 1. Using non-phylogenetic methods: Dice-Sorensen, Bray-Curtis, Morisita-Horn
## 2. Using Unifrac: weighted and unweighted

## 1: export Physeq object and calculating non-phylogenetic distances in Vegan
nonPhyloDists.function <- function(Phydata, BetaTest, BinaryVal, BetaMethod, Pid) {
  newPhy <- subset_samples(Phydata, SampleType != 'control')   ## drops Control samples from analysis
  data_mat <- as(otu_table(newPhy), "matrix")
  tmp.meta <- data.frame(tinymeta) %>% filter(SampleID %in% row.names(data_mat))
  row.names(tmp.meta) <- tmp.meta$SampleID
  betadist <- vegdist(x = data_mat, method = BetaTest, binary = BinaryVal)
  tmp.adonis <- adonis(betadist ~ CollectionMonth * Site * Roost, data = tmp.meta)
  tmp.adonis <- as.data.frame(tmp.adonis$aov.tab)
  tmp.out <- data.frame(tmp.adonis, BetaMethod, Pid)
  tmp.out$TestGroup <- row.names(tmp.out)
  tmp.out
}

dice.p100.adonis <- nonPhyloDists.function(rphy_p100, "bray", TRUE, "dice", "pid=100")
bray.p100.adonis <- nonPhyloDists.function(rphy_p100, "bray", FALSE, "bray", "pid=100")
mori.p100.adonis <- nonPhyloDists.function(rphy_p100, "morisita", FALSE, "morisita", "pid=100")

dice.p99.adonis <- nonPhyloDists.function(rphy_p99, "bray", TRUE, "dice", "pid=99")
bray.p99.adonis <- nonPhyloDists.function(rphy_p99, "bray", FALSE, "bray", "pid=99")
mori.p99.adonis <- nonPhyloDists.function(rphy_p99, "morisita", FALSE, "morisita", "pid=99")

dice.p98.adonis <- nonPhyloDists.function(rphy_p98, "bray", TRUE, "dice", "pid=98")
bray.p98.adonis <- nonPhyloDists.function(rphy_p98, "bray", FALSE, "bray", "pid=98")
mori.p98.adonis <- nonPhyloDists.function(rphy_p98, "morisita", FALSE, "morisita", "pid=98")

nonPhylo_adonis_all <- rbind(dice.p100.adonis, bray.p100.adonis, mori.p100.adonis, dice.p99.adonis, bray.p99.adonis, mori.p99.adonis, dice.p98.adonis, bray.p98.adonis, mori.p98.adonis)
write.csv(nonPhylo_adonis_all, file = "~/Repos/mysosoup/data/adonis.nonPhylo.csv", quote = FALSE, row.names = FALSE)
rm(dice.p100.adonis, bray.p100.adonis, mori.p100.adonis, dice.p99.adonis, bray.p99.adonis, mori.p99.adonis, dice.p98.adonis, bray.p98.adonis, mori.p98.adonis)


## 2: import rooted tree; drop outliers; calculate phylogenetic-informed distances in Phyloseq; export distances to run Adonis
## import tree object
pid100_tree <- read.tree(file = "~/Repos/mysosoup/data/trees/ASV_famFilt.nwk")

## 'BadSamples' dropped based from dataset for weighted Unifract to work
## These samples were identified as outliers from the `physeqNMDS_unifrac.R` script
'%!in%' <- function(x,y)!('%in%'(x,y))
BadSamples <- c('ExtractionNTC4S28', '7272017EGB9', '9152017EGC3', '6212017HBB4')

## function to drop group of samples
nonPhyloDists.function <- function(Phydata, Pid) {
  rphy_wTree <- merge_phyloseq(Phydata, pid100_tree)
  newPhy <- subset_samples(rphy_wTree, SampleType != 'control')   ## drops Control samples from analysis
  newPhy <- subset_samples(newPhy, SampleID %!in% BadSamples)
  rm(rphy_wTree)
  tmp.meta <- data.frame(tinymeta) %>% filter(SampleID %in% sample_names(newPhy))
  uUF_dist <- UniFrac(newPhy, weighted = FALSE)
  u.adonis <- adonis(uUF_dist ~ CollectionMonth * Site * Roost, data = tmp.meta)
  u.adonis <- as.data.frame(u.adonis$aov.tab)
  BetaMethod = 'Unifrac'
  Weight = TRUE
  u.out <- data.frame(u.adonis, BetaMethod, Weight, Pid)
  u.out$TestGroup <- row.names(u.out)
  wUF_dist <- UniFrac(newPhy, weighted = TRUE)
  w.adonis <- adonis(wUF_dist ~ CollectionMonth * Site * Roost, data = tmp.meta)
  w.adonis <- as.data.frame(w.adonis$aov.tab)
  BetaMethod = 'Unifrac'
  Weight = FALSE
  w.out <- data.frame(w.adonis, BetaMethod, Weight, Pid)
  w.out$TestGroup <- row.names(w.out)
  rbind(u.out, w.out)
}
  
uni_p100_df <- nonPhyloDists.function(rphy_p100, "pid=100")
uni_p99_df <- nonPhyloDists.function(rphy_p100, "pid=99")
uni_p98_df <- nonPhyloDists.function(rphy_p100, "pid=98")

all_uni_adonis <- rbind(uni_p100_df, uni_p99_df, uni_p98_df)
write.csv(all_uni_adonis, file="~/Repos/mysosoup/data/adonis.Phylo.csv", quote = FALSE, row.names = FALSE)
