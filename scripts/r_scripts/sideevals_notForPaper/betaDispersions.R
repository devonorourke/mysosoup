library(tidyverse)
library(phyloseq)
library(reshape2)
library(vegan)
library(ape)
library(qiime2R)

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

## arbitrarily picking 9000 because it appears to capture most data without going into left tail of distribution
rphy_p100 = rarefy_even_depth(phy_p100, sample.size=9000, replace = FALSE, rngseed = 1000)
rphy_p99 = rarefy_even_depth(phy_p99, sample.size=9000, replace = FALSE, rngseed = 1000)
rphy_p98 = rarefy_even_depth(phy_p98, sample.size=9000, replace = FALSE, rngseed = 1000)

rm(phy_p100, phy_p99, phy_p98)

## set up 2 functions to calculate distances:
## 1. Using non-phylogenetic methods: Dice-Sorensen, Bray-Curtis, Morisita-Horn
## 2. Using Unifrac: weighted and unweighted

## 1: export Physeq object and calculating non-phylogenetic distances in Vegan
nonPhyloDisp.function <- function(Phydata, BetaTest, BinaryVal, DistMethod, Pid) {
  #newPhy <- subset_samples(rphy_p100, SampleType != 'control')   ## drops Control samples from analysis
  newPhy <- subset_samples(Phydata, SampleType != 'control')   ## drops Control samples from analysis
  data_mat <- as(otu_table(newPhy), "matrix")
  tmp.meta <- data.frame(tinymeta) %>% filter(SampleID %in% row.names(data_mat))
  row.names(tmp.meta) <- tmp.meta$SampleID
  tmp.meta$Group <- paste(tmp.meta$Labeler, tmp.meta$CollectionMonth, sep = "-")
  betadist <- vegdist(x = data_mat, method = BetaTest, binary = BinaryVal)
  groups <- as.character(tmp.meta$Group)
  tmp.disp <- betadisper(betadist, groups, type = c("median"))
  rm(betadist)
  tmp.disp_sco <- as.data.frame(scores(tmp.disp, display = "sites"))
  tmp.disp_sco$bigID <- row.names(tmp.disp_sco)
  row.names(tmp.disp_sco) <- NULL
  colnames(tmp.disp_sco) <- c("pc1end", "pc2end", "SampleID")
  tmp.disp_sco <- merge(tmp.disp_sco, tinymeta)
  tmp.disp_sco$Group <- paste(tmp.disp_sco$Site, tmp.disp_sco$Roost, tmp.disp_sco$CollectionMonth, sep = "-")
  tmp.disp_centroids <- as.data.frame(scores(tmp.disp, display="centroids"))
  tmp.disp_centroids$SeqID <- row.names(tmp.disp_centroids)
  row.names(tmp.disp_centroids) <- NULL
  colnames(tmp.disp_centroids) <- c("pc1start", "pc2start", "Group")
  tmp.disp.plot <- merge(tmp.disp_sco, tmp.disp_centroids)
  tmp.disp.plot <- tmp.disp.plot %>% mutate(DistMethod=DistMethod) %>% mutate(Pid)
  tmp.disp.plot
}

dice.p100.disp <- nonPhyloDisp.function(rphy_p100, "bray", TRUE, "dice", "pid=100")
bray.p100.disp <- nonPhyloDisp.function(rphy_p100, "bray", FALSE, "bray", "pid=100")
mori.p100.disp <- nonPhyloDisp.function(rphy_p100, "morisita", FALSE, "morisita", "pid=100")

dice.p99.disp <- nonPhyloDisp.function(rphy_p99, "bray", TRUE, "dice", "pid=99")
bray.p99.disp <- nonPhyloDisp.function(rphy_p99, "bray", FALSE, "bray", "pid=99")
mori.p99.disp <- nonPhyloDisp.function(rphy_p99, "morisita", FALSE, "morisita", "pid=99")

dice.p98.disp <- nonPhyloDisp.function(rphy_p98, "bray", TRUE, "dice", "pid=98")
bray.p98.disp <- nonPhyloDisp.function(rphy_p98, "bray", FALSE, "bray", "pid=98")
mori.p98.disp <- nonPhyloDisp.function(rphy_p98, "morisita", FALSE, "morisita", "pid=98")

nonPhylo_disp_all <- rbind(dice.p100.disp, bray.p100.disp, mori.p100.disp, dice.p99.disp, bray.p99.disp, mori.p99.disp, dice.p98.disp, bray.p98.disp, mori.p98.disp)
write.csv(nonPhylo_disp_all, file = "~/Repos/mysosoup/data/disp.nonPhylo.csv", quote = FALSE, row.names = FALSE)
rm(dice.p100.disp, bray.p100.disp, mori.p100.disp, dice.p99.disp, bray.p99.disp, mori.p99.disp, dice.p98.disp, bray.p98.disp, mori.p98.disp)

## set levels for plot
nonPhylo_disp_all$DistMethod <- factor(nonPhylo_disp_all$DistMethod, levels = c("dice", "bray", "morisita"))
nonPhylo_disp_all$CollectionMonth <- factor(nonPhylo_disp_all$CollectionMonth, levels = c("June", "July", "September"))
nonPhylo_disp_all$Pid <- factor(nonPhylo_disp_all$Pid, levels = c("pid=100", "pid=99", "pid=98"))

## colors for plot
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)

## plot dispersions; points are samples
## save as nonPhylo_Dispersions; export at 1200x1200
ggplot() + 
  geom_point(data=nonPhylo_disp_all, aes(x=pc1end, y=pc2end, color=CollectionMonth), size=2, alpha=0.8) +
  scale_color_manual(values=v3pal) + 
  geom_segment(data=nonPhylo_disp_all, aes(x=pc1end, y=pc2end, xend=pc1start, yend=pc2start), alpha=0.12) +
  facet_grid(Pid ~ DistMethod) +
  labs(title="", color="", x="PCoA1", y="PCoA2", shape="", color="") +
  theme_devon() +
  theme(legend.position = "top")


## 2: import rooted tree; drop outliers; calculate phylogenetic-informed distances in Phyloseq; export distances to run Adonis
## import tree object
pid100_tree <- read.tree(file = "~/Repos/mysosoup/data/trees/ASV_famFilt.nwk")

## 'BadSamples' dropped based from dataset for weighted Unifract to work
## These samples were identified as outliers from the `physeqNMDS_unifrac.R` script
'%!in%' <- function(x,y)!('%in%'(x,y))
BadSamples <- c('ExtractionNTC4S28', '7272017EGB9', '9152017EGC3', '6212017HBB4')

## function to drop group of samples
PhyloDisp.function <- function(Phydata, BinaryVal, DistMethod, Pid) {
  rphy_wTree <- merge_phyloseq(Phydata, pid100_tree)
  newPhy <- subset_samples(rphy_wTree, SampleType != 'control')   ## drops Control samples from analysis
  newPhy <- subset_samples(newPhy, SampleID %!in% BadSamples)
  rm(rphy_wTree)
  tmp.meta <- data.frame(tinymeta) %>% filter(SampleID %in% sample_names(newPhy))
  tmp.meta$Group <- paste(tmp.meta$Labeler, tmp.meta$CollectionMonth, sep = "-")
  UF_dist <- UniFrac(newPhy, weighted = BinaryVal)
  groups <- as.character(tmp.meta$Group)
  tmp.disp <- betadisper(UF_dist, groups, type = c("median"))
  rm(UF_dist)
  tmp.disp_sco <- as.data.frame(scores(tmp.disp, display = "sites"))
  tmp.disp_sco$bigID <- row.names(tmp.disp_sco)
  row.names(tmp.disp_sco) <- NULL
  colnames(tmp.disp_sco) <- c("pc1end", "pc2end", "SampleID")
  tmp.disp_sco <- merge(tmp.disp_sco, tinymeta)
  tmp.disp_sco$Group <- paste(tmp.disp_sco$Site, tmp.disp_sco$Roost, tmp.disp_sco$CollectionMonth, sep = "-")
  tmp.disp_centroids <- as.data.frame(scores(tmp.disp, display="centroids"))
  tmp.disp_centroids$SeqID <- row.names(tmp.disp_centroids)
  row.names(tmp.disp_centroids) <- NULL
  colnames(tmp.disp_centroids) <- c("pc1start", "pc2start", "Group")
  tmp.disp.plot <- merge(tmp.disp_sco, tmp.disp_centroids)
  tmp.disp.plot <- tmp.disp.plot %>% mutate(DistMethod=DistMethod) %>% mutate(Pid=Pid)
  tmp.disp.plot
}
  
wUF.p100.disp <- PhyloDisp.function(rphy_p100, TRUE, "wUnifrac", "pid=100")
wUF.p99.disp <- PhyloDisp.function(rphy_p99, TRUE, "wUnifrac", "pid=99")
wUF.p98.disp <- PhyloDisp.function(rphy_p98, TRUE, "wUnifrac", "pid=98")
uUF.p100.disp <- PhyloDisp.function(rphy_p100, FALSE, "uUnifrac", "pid=100")
uUF.p99.disp <- PhyloDisp.function(rphy_p99, FALSE, "uUnifrac", "pid=99")
uUF.p98.disp <- PhyloDisp.function(rphy_p98, FALSE, "uUnifrac", "pid=98")

UF_all.disp <- rbind(wUF.p100.disp, wUF.p99.disp, wUF.p98.disp, uUF.p100.disp, uUF.p99.disp, uUF.p98.disp)
rm(wUF.p100.disp, wUF.p99.disp, wUF.p98.disp, uUF.p100.disp, uUF.p99.disp, uUF.p98.disp)


## set levels for plot
UF_all.disp$DistMethod <- factor(UF_all.disp$DistMethod, levels = c("uUnifrac", "wUnifrac"))
UF_all.disp$CollectionMonth <- factor(UF_all.disp$CollectionMonth, levels = c("June", "July", "September"))
UF_all.disp$Pid <- factor(UF_all.disp$Pid, levels = c("pid=100", "pid=99", "pid=98"))

## save as Phylo_Dispersions; export at 900x900
ggplot() + 
  geom_point(data=UF_all.disp, aes(x=pc1end, y=pc2end, color=CollectionMonth), size=2, alpha=0.8) +
  scale_color_manual(values=v3pal) + 
  geom_segment(data=UF_all.disp, aes(x=pc1end, y=pc2end, xend=pc1start, yend=pc2start), alpha=0.12) +
  facet_grid(Pid ~ DistMethod) +
  labs(title="", color="", x="PCoA1", y="PCoA2", shape="", color="") +
  theme_devon() +
  theme(legend.position = "top")


## combine all dispersion data for a single plot:
## save as all_Dispersions; export at 1200-900
all.disp <- rbind(UF_all.disp, nonPhylo_disp_all)
ggplot() + 
  geom_point(data=all.disp, aes(x=pc1end, y=pc2end, color=CollectionMonth), size=2, alpha=0.8) +
  scale_color_manual(values=v3pal) + 
  geom_segment(data=all.disp, aes(x=pc1end, y=pc2end, xend=pc1start, yend=pc2start), alpha=0.12) +
  facet_grid(Pid ~ DistMethod) +
  labs(title="", color="", x="PCoA1", y="PCoA2", shape="", color="") +
  theme_devon() +
  theme(legend.position = "top")


#### ----- Stats for Beta Dispersions

PhyloDispStats.function <- function(Phydata, BinaryVal, DistMethod, Pid) {
  #rphy_wTree <- merge_phyloseq(rphy_p100, pid100_tree)
  rphy_wTree <- merge_phyloseq(Phydata, pid100_tree)
  newPhy <- subset_samples(rphy_wTree, SampleType != 'control')   ## drops Control samples from analysis
  newPhy <- subset_samples(newPhy, SampleID %!in% BadSamples)
  rm(rphy_wTree)
  tmp.meta <- data.frame(tinymeta) %>% filter(SampleID %in% sample_names(newPhy))
  tmp.meta$Group <- paste(tmp.meta$Labeler, tmp.meta$CollectionMonth, sep = "-")
  UF_dist <- UniFrac(newPhy, weighted = BinaryVal)
  groups <- as.character(tmp.meta$Group)
  tmp.disp <- betadisper(UF_dist, groups, type = c("median"))
  aov_out <- data.frame(anova(tmp.disp))
  bdips.dist <- data.frame(tmp.disp$distances, tmp.disp$group)
  colnames(bdips.dist) <- c("distances", "group")
  bdips.dist$sample <- row.names(bdips.dist)
  row.names(bdips.dist) <- NULL
  bdips.dist <- bdips.dist %>% separate(., col = "group", into = c("Site", "Roost", "CollectionMonth"), sep = "-")
  bdips.anova <- aov(distances ~ Site * Roost * CollectionMonth, data = bdips.dist)
  aov_out2 <- summary(bdips.anova)
  list(aov_out, aov_out2)
}

## generate anova calculations treating entire unit as one factor, or as separate factors
## access either anova with: 
## uUF.p100.dispStats.list[1] ... or ...
## uUF.p100.dispStats.list[2]

uUF.p100.dispStats.list <- PhyloDispStats.function(rphy_p100, FALSE, 'uUnifrac', 'pid=100')
uUF.p99.dispStats.list <- PhyloDispStats.function(rphy_p99, FALSE, 'uUnifrac', 'pid=99')
uUF.p98.dispStats.list <- PhyloDispStats.function(rphy_p98, FALSE, 'uUnifrac', 'pid=98')
wUF.p100.dispStats.list <- PhyloDispStats.function(rphy_p100, TRUE, 'wUnifrac', 'pid=100')
wUF.p99.dispStats.list <- PhyloDispStats.function(rphy_p99, TRUE, 'wUnifrac', 'pid=99')
wUF.p98.dispStats.list <- PhyloDispStats.function(rphy_p98, TRUE, 'wUnifrac', 'pid=98')

out1a <- data.frame(uUF.p100.dispStats.list[1])
out1b <- data.frame(wUF.p100.dispStats.list[1])
out1c <- data.frame(uUF.p99.dispStats.list[1])
out1d <- data.frame(wUF.p99.dispStats.list[1])
out1e <- data.frame(uUF.p98.dispStats.list[1])
out1f <- data.frame(wUF.p98.dispStats.list[1])

write.csv(out1a, file="~/Repos/mysosoup/data/text_tables/uUF-p100.dispAnova.single.csv", quote = FALSE, row.names = TRUE)
write.csv(out1b, file="~/Repos/mysosoup/data/text_tables/wUF-p100.dispAnova.single.csv", quote = FALSE, row.names = TRUE)
write.csv(out1c, file="~/Repos/mysosoup/data/text_tables/uUF-p99.dispAnova.single.csv", quote = FALSE, row.names = TRUE)
write.csv(out1d, file="~/Repos/mysosoup/data/text_tables/wUF-p99.dispAnova.single.csv", quote = FALSE, row.names = TRUE)
write.csv(out1e, file="~/Repos/mysosoup/data/text_tables/uUF-p98.dispAnova.single.csv", quote = FALSE, row.names = TRUE)
write.csv(out1f, file="~/Repos/mysosoup/data/text_tables/wUF-p98.dispAnova.single.csv", quote = FALSE, row.names = TRUE)

out2a <- do.call(rbind.data.frame, uUF.p100.dispStats.list[[2]][1])
out2b <- do.call(rbind.data.frame, uUF.p99.dispStats.list[[2]][1])
out2c <- do.call(rbind.data.frame, uUF.p98.dispStats.list[[2]][1])
out2d <- do.call(rbind.data.frame, wUF.p100.dispStats.list[[2]][1])
out2e <- do.call(rbind.data.frame, wUF.p99.dispStats.list[[2]][1])
out2f <- do.call(rbind.data.frame, wUF.p98.dispStats.list[[2]][1])

write.csv(out2a, file="~/Repos/mysosoup/data/text_tables/uUF-p100.dispAnova.multiple.csv", quote = FALSE, row.names = TRUE)
write.csv(out2b, file="~/Repos/mysosoup/data/text_tables/uUF-p99.dispAnova.multiple.csv", quote = FALSE, row.names = TRUE)
write.csv(out2c, file="~/Repos/mysosoup/data/text_tables/uUF-p98.dispAnova.multiple.csv", quote = FALSE, row.names = TRUE)
write.csv(out2d, file="~/Repos/mysosoup/data/text_tables/wUF-p100.dispAnova.multiple.csv", quote = FALSE, row.names = TRUE)
write.csv(out2e, file="~/Repos/mysosoup/data/text_tables/wUF-p99.dispAnova.multiple.csv", quote = FALSE, row.names = TRUE)
write.csv(out2f, file="~/Repos/mysosoup/data/text_tables/wUF-p98.dispAnova.multiple.csv", quote = FALSE, row.names = TRUE)