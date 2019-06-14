library(tidyverse)
library(qiime2R)
library(reshape2)
library(phyloseq)
library(ape)

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
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE)
tinymeta <- meta %>% select(SampleID, Roost, CollectionMonth, SampleType, Site, ContamArea)
tinymeta$Site <- ifelse(tinymeta$Site == "Egner", gsub("Egner", "EN", tinymeta$Site), tinymeta$Site)
tinymeta$Site <- ifelse(tinymeta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", tinymeta$Site), tinymeta$Site)
tinymeta$CollectionMonth[is.na(tinymeta$CollectionMonth)] <- "control"
tinymeta$SiteMonth <- paste(tinymeta$Site, tinymeta$CollectionMonth, sep="-")
tinymeta$SiteMonth <- ifelse(tinymeta$SiteMonth == "control-control", gsub("control-control", "control", tinymeta$SiteMonth), tinymeta$SiteMonth)
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "6", gsub("6", "June", tinymeta$CollectionMonth), tinymeta$CollectionMonth)
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "7", gsub("7", "July", tinymeta$CollectionMonth), tinymeta$CollectionMonth)
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "9", gsub("9", "September", tinymeta$CollectionMonth), tinymeta$CollectionMonth)

## import taxonomy info
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"


## create physeq metadata and taxonomy data for analyses:
row.names(meta) <- meta$SampleID
phy_meta <- sample_data(meta)
row.names(taxa) <- taxa$ASVid
phy_tax <- as.matrix(taxa) %>% tax_table(.)

## import sequence data and combine Taxa, Meta, and read data into phyloseq object
## performed for each of 2 datasets (w and w/out ASVs in NTCs)
physeqImport <- function(qzaPath){
  featuretable <- read_qza(qzaPath)
  mat.tmp <- featuretable$data
  OTU <- otu_table(mat.tmp, taxa_are_rows = TRUE) ## import as physeq object 
  phyloseq(OTU, phy_tax, phy_meta)
}

## filepaths to import QZA files from:
noNTCasv_path <- "~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.noNTCasvs-filt.rarefied-table.qza"
wNTCasv_path <- "~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza"

## applying Physeq import function:
phy_wNTCasv <- physeqImport(wNTCasv_path)
phy_noNTCasv <- physeqImport(noNTCasv_path)

#rm(phy_tax, phy_meta)
nsamples(phy_wNTCasv)  ## 288 samples remain (from 304, 11 of which were NTCs)
nsamples(phy_noNTCasv)  ## 265 samples remain (includes some NTC samples)
ntaxa(phy_wNTCasv) ## 2583 taxa remain
ntaxa(phy_noNTCasv) ## 2379 taxa remain (fewer, in part, because we dropped ASVs from NTCs)

## Beta diversity with phylogeny
## import tree object
tree <- read.tree(file = "~/Repos/mysosoup/data/trees/rooted-tree.wNTCasvs.nwk")

## calculate distances with unifrac
BetaPhyFunction <- function(Phydata, FiltTable) {
  rphy_wTree <- merge_phyloseq (Phydata, tree)
  wUF_ord = ordinate(rphy_wTree, method="NMDS", distance="unifrac", weighted=TRUE)
  wUF_df <- data.frame(wUF_ord$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="wu") %>% mutate(FiltTable=FiltTable)
  rm(wUF_ord)
  uUF_ord = ordinate(rphy_wTree, method="NMDS", distance="unifrac", weighted=FALSE)
  uUF_df <- data.frame(uUF_ord$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="uu") %>% mutate(FiltTable=FiltTable)
  rm(uUF_ord)
  out <- rbind(wUF_df, uUF_df)
  out2 <- merge(tinymeta, out)
  out2
}

wNTCasv_uni_df <- BetaPhyFunction(phy_wNTCasv, "wNTCasv")
  ## high stress (~ 0.19); no solution reached
noNTCasv_uni_df <- BetaPhyFunction(phy_noNTCasv, "noNTCasv")
  ## high stress (~ 0.18); solution reached
ALL_uni_df <- rbind(wNTCasv_uni_df, noNTCasv_uni_df)


## plot; save as 'nmds.unifrac_allTables'; export at 800x800
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)
ALL_uni_df$CollectionMonth <- factor(ALL_uni_df$CollectionMonth, levels=c("June", "July", "September", "control"))
ggplot(ALL_uni_df, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site, label=SampleID)) +
  geom_point() +
  scale_color_manual(values=c(v3pal, "blue")) +
  facet_grid(FiltTable ~ Measure) +
  theme_devon() +
  theme(legend.position = "top")

## individual plots; save as 'nmds.unifrac_wNTCasv_wControlsShown'; export at 800x475
wNTCasv_uni_df$CollectionMonth <- factor(wNTCasv_uni_df$CollectionMonth, levels=c("June", "July", "September", "control"))
ggplot(wNTCasv_uni_df, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site)) +
  geom_point() +
  scale_color_manual(values=c(v3pal, "blue")) +
  facet_grid(FiltTable ~ Measure) +
  theme_devon() +
  theme(legend.position = "top")

## individual plots; save as 'nmds.unifrac_noNTCasv'; export at 800x475
noNTCasv_uni_df$CollectionMonth <- factor(noNTCasv_uni_df$CollectionMonth, levels=c("June", "July", "September"))
ggplot(noNTCasv_uni_df, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site)) +
  geom_point() +
  scale_color_manual(values=c(v3pal, "blue")) +
  facet_grid(FiltTable ~ Measure) +
  theme_devon() +
  theme(legend.position = "top")


#### repeat plots, but calculate distances with non-phylogenetic fucntions 
## going to find weird outliers here -- need to remove samples
## executing functions independently to calculate distances so we can remove odd samples unique to each dataset

tmp_wTree <- merge_phyloseq (phy_wNTCasv, tree)
tmp_ord = ordinate(tmp_wTree, method="NMDS", distance="bray", binary=TRUE, transformation=FALSE)
ds_plot <- data.frame(tmp_ord$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="ds") %>% mutate(FiltTable="wNTCasv")


## consider dropping that one outlier?
'%!in%' <- function(x,y)!('%in%'(x,y))
wNTC_BadSamples <- c('7272017EGC1', '7272017HBA6', 'ExtractionNTC11S115')

nonPhyBetaFunction_wNTC <- function(Phydata, FiltTable) {
  rphy_wTree <- merge_phyloseq (Phydata, tree)
  rphy_wTree <- subset_samples(rphy_wTree, SampleID %!in% wNTC_BadSamples)
  dist_bc <- phyloseq::distance(rphy_wTree, "bray")
  nmds_bc <- ordinate(rphy_wTree, method = "NMDS", distance = dist_bc)
  dat_bc <- data.frame(nmds_bc$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="bc") %>% mutate(FiltTable=FiltTable)
  dist_ds <- phyloseq::distance(rphy_wTree, "bray", binary = TRUE)
  nmds_ds <- ordinate(rphy_wTree, method = "NMDS", distance = dist_ds)
  dat_ds <- data.frame(nmds_ds$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="ds") %>% mutate(FiltTable=FiltTable)
  dist_mh <- phyloseq::distance(rphy_wTree, "morisita", binary = FALSE)
  nmds_mh <- ordinate(rphy_wTree, method = "NMDS", distance = dist_mh)
  dat_mh <- data.frame(nmds_mh$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="mh") %>% mutate(FiltTable=FiltTable)
  out <- rbind(dat_bc, dat_ds, dat_mh)
  out <- merge(out, tinymeta)
  out
}

wNTCasv_plot <- nonPhyBetaFunction_wNTC(phy_wNTCasv, "wNTCasv")
wNTCasv_plot$ContamArea[is.na(wNTCasv_plot$ContamArea)] <- "isolate"

## wNTC plot
wNTCasv_plot$CollectionMonth <- factor(wNTCasv_plot$CollectionMonth, levels=c("June", "July", "September", "control"))
wNTCasv_plot$Measure <- factor(wNTCasv_plot$Measure, levels=c("ds", "bc", "mh"))

## plot with control samples highlighted; save as 'wNTCasvs_NMDS_wControls_byMonth' export at 1000x475
ggplot(wNTCasv_plot, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site)) +
  geom_point(data = wNTCasv_plot %>% filter(SampleType == "sample"), aes(color=CollectionMonth), alpha=0.7) +
  geom_point(data = wNTCasv_plot %>% filter(SampleType == "control"), color="blue", size=2) +
  scale_color_manual(values=c(v3pal, "gray40")) +
  facet_grid( ~ Measure) +
  theme_devon() +
  theme(legend.position = "top")

## same plot, but coloring by contamination area; save as 'wNTCasvs_NMDS_wControls_byContamArea'; export at 1000x475
ggplot(wNTCasv_plot, aes(x=MDS1, y=MDS2, color=Site, shape=ContamArea)) +
  geom_point(data = wNTCasv_plot %>% filter(ContamArea != TRUE)) +
  geom_point(data = wNTCasv_plot %>% filter(ContamArea == TRUE), size=2.5) +
  scale_color_manual(values=c("red", "navy", "gray40")) +
  facet_grid( ~ Measure) +
  theme_devon() +
  theme(legend.position = "top")


## noNTC next ...
noNTC_BadSamples <- c('7272017EGC1', '7272017EGC2', '7272017HBA10', '7272017HBD10', '7272017EGB9', '7272017HBA6', '7272017HBA6', '6212017EGD8')

nonPhyBetaFunction_noNTC <- function(Phydata, FiltTable) {
  rphy_wTree <- merge_phyloseq (Phydata, tree)
  rphy_wTree <- subset_samples(rphy_wTree, SampleID %!in% noNTC_BadSamples)
  dist_bc <- phyloseq::distance(rphy_wTree, "bray")
  nmds_bc <- ordinate(rphy_wTree, method = "NMDS", distance = dist_bc)
  dat_bc <- data.frame(nmds_bc$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="bc") %>% mutate(FiltTable=FiltTable)
  dist_ds <- phyloseq::distance(rphy_wTree, "bray", binary = TRUE)
  nmds_ds <- ordinate(rphy_wTree, method = "NMDS", distance = dist_ds)
  dat_ds <- data.frame(nmds_ds$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="ds") %>% mutate(FiltTable=FiltTable)
  dist_mh <- phyloseq::distance(rphy_wTree, "morisita", binary = FALSE)
  nmds_mh <- ordinate(rphy_wTree, method = "NMDS", distance = dist_mh)
  dat_mh <- data.frame(nmds_mh$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="mh") %>% mutate(FiltTable=FiltTable)
  out <- rbind(dat_bc, dat_ds, dat_mh)
  out <- merge(out, tinymeta)
  out
}

noNTCasv_plot <- nonPhyBetaFunction_noNTC(phy_noNTCasv, "noNTCasv")

## noNTC plot
noNTCasv_plot$CollectionMonth <- factor(noNTCasv_plot$CollectionMonth, levels=c("June", "July", "September"))
noNTCasv_plot$Measure <- factor(noNTCasv_plot$Measure, levels = c("ds", "bc", "mh"))

## save as 'noNTCasvs_NMDS_bySiteMonth'; export at 1000x475
ggplot(noNTCasv_plot, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site, label=SampleID)) +
  geom_point() +
  scale_color_manual(values=c(v3pal, "gray40")) +
  facet_grid( ~ Measure) +
  theme_devon() +
  theme(legend.position = "top")


## These data suggest that the ASVs we're filtering aren't associated with Site nor Month nor well location (ContamArea)
## We're going to keep all the ASVs for the final beta diversity plot, but remove the contaminated samples first:
## consider dropping that one outlier?

wNTC_dropSamples <- c('blankS39', 'ExtractionNTC11S115', "ExtractionNTC1S1", "ExtractionNTC2S10",
                     "ExtractionNTC3S19", "ExtractionNTC4S28", "ExtractionNTC7S55", "ExtractionNTC8S64",
                     '7272017EGC1', '7272017HBA6')

nonPhyBetaFunction_wNTC_NTCdropd <- function(Phydata, FiltTable) {
  rphy_wTree <- merge_phyloseq (Phydata, tree)
  rphy_wTree <- subset_samples(rphy_wTree, SampleID %!in% wNTC_dropSamples)
  dist_bc <- phyloseq::distance(rphy_wTree, "bray")
  nmds_bc <- ordinate(rphy_wTree, method = "NMDS", distance = dist_bc)
  dat_bc <- data.frame(nmds_bc$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="bc") %>% mutate(FiltTable=FiltTable)
  dist_ds <- phyloseq::distance(rphy_wTree, "bray", binary = TRUE)
  nmds_ds <- ordinate(rphy_wTree, method = "NMDS", distance = dist_ds)
  dat_ds <- data.frame(nmds_ds$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="ds") %>% mutate(FiltTable=FiltTable)
  dist_mh <- phyloseq::distance(rphy_wTree, "morisita", binary = FALSE)
  nmds_mh <- ordinate(rphy_wTree, method = "NMDS", distance = dist_mh)
  dat_mh <- data.frame(nmds_mh$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="mh") %>% mutate(FiltTable=FiltTable)
  wUF_ord = ordinate(rphy_wTree, method="NMDS", distance="unifrac", weighted=TRUE)
  wUF_df <- data.frame(wUF_ord$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="wu") %>% mutate(FiltTable=FiltTable)
  rm(wUF_ord)
  uUF_ord = ordinate(rphy_wTree, method="NMDS", distance="unifrac", weighted=FALSE)
  uUF_df <- data.frame(uUF_ord$points) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="uu") %>% mutate(FiltTable=FiltTable)
  rm(uUF_ord)
  out <- rbind(dat_bc, dat_ds, dat_mh, wUF_df, uUF_df)
  out <- merge(out, tinymeta)
  out
}


wNTCasv_noNTC_plot <- nonPhyBetaFunction_wNTC_NTCdropd(phy_wNTCasv, "wNTCasv")
wNTCasv_noNTC_plot$CollectionMonth <- factor(wNTCasv_noNTC_plot$CollectionMonth, levels=c("June", "July", "September"))
wNTCasv_noNTC_plot$Measure <- factor(wNTCasv_noNTC_plot$Measure, levels=c("ds", "bc", "mh", "uu", "wu"))

## plot with control samples highlighted; save as 'wNTCasvs_NMDS_noControls_byMonth' export at 1000x475
ggplot(wNTCasv_noNTC_plot, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site)) +
  geom_point(data = wNTCasv_noNTC_plot %>% filter(SampleType == "sample"), aes(color=CollectionMonth)) +
  scale_color_manual(values=c(v3pal)) +
  facet_wrap( ~ Measure, nrow = 2) +
  theme_devon() +
  theme(legend.position = "top")

