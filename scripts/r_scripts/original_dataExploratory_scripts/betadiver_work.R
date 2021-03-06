## see here for more helpful physeq tricks: https://raw.githubusercontent.com/MadsAlbertsen/ampvis/gh-pages/examples/ampvis_guide.Rmd
## more tricks: http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html#multidimensional-scaling-mds-aka-pcoa

library(tidyverse)
library(qiime2R)
library(reshape2)
library(phyloseq)
library(ape)
library(vegan)

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

########################################################################
## Part 1 == data wrangling: import meta, taxa, read data; 
##        == create phyloseq object; rarefy object
########################################################################

## import metadata 
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE)
meta <- meta %>% 
  select(SampleID, Roost, CollectionMonth, SampleType, Site, ContamArea) %>% 
  rename(Month = CollectionMonth)
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$SiteMonth <- paste(meta$Site, meta$Month, sep="-")
meta$Month <- ifelse(meta$Month == "6", gsub("6", "June", meta$Month), meta$Month)
meta$Month <- ifelse(meta$Month == "7", gsub("7", "July", meta$Month), meta$Month)
meta$Month <- ifelse(meta$Month == "9", gsub("9", "September", meta$Month), meta$Month)
meta <- meta %>% filter(SampleType == "sample")

## import taxonomy info
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"
row.names(taxa) <- taxa$ASVid


## create physeq metadata and taxonomy data for analyses:
row.names(meta) <- meta$SampleID
phy_meta <- sample_data(meta)
row.names(taxa) <- taxa$ASVid
phy_tax <- as.matrix(taxa) %>% tax_table(.)

## import sequence data from QZA artifact, filter out reads all non-Arthropod reads without (at least) Family-rank info, .. 
## ..filter out remaining ASVs present exclusively in control samples, then import as physeq object
## download qza file:
# download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza", "tmp.qza")
# then set PATH to whatever directory needed:
# qzapath = "PATH/TO/tmp.qza"

## alternatively run from local: qzapath = "~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza"
featuretable <- read_qza(qzapath)
mat.tmp <- featuretable$data
rm(featuretable)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
rm(df.tmp) 
colnames(tmp) <- c("ASVid", "SampleID", "Reads")
df.tmp <- merge(tmp, taxa)
rm(tmp)
df.tmp <- merge(df.tmp, meta)
mat.tmp <- dcast(data = df.tmp, formula = SampleID ~ ASVid, value.var='Reads', fill = 0)
row.names(mat.tmp) <- mat.tmp$SampleID
mat.tmp$SampleID <- NULL
OTU <- otu_table(mat.tmp, taxa_are_rows = FALSE) ## import as physeq object 
phydat <- phyloseq(OTU, phy_tax, phy_meta)
rm(mat.tmp, phy_tax, OTU, qzapath, phy_meta)

## sanity check:
nsamples(phydat)  ## 277 samples remain (from 304, 11 of which were NTCs)
ntaxa(phydat) ## 2573 taxa remain
sample_names(phydat)  ## all contaminant samples removed
  ## all good.

########################################################################
## Part 2 == calculate distances using 5 unique distance measures
##        == ordinate and plot with 95% confidence intervals ~ Site * Month
########################################################################

## import tree object to for phylogenetic distance measures
## download qza file:
# download.file("https://github.com/devonorourke/mysosoup/raw/master/data/trees/rooted-tree.wNTCasvs.nwk", "tree.nwk")
# then set PATH to whatever directory needed:
# tree = "PATH/TO/tree.nwk"

## alternatively run from local: tree <- read.tree(file = "~/Repos/mysosoup/data/trees/rooted-tree.wNTCasvs.nwk")

## add tree info to physeq object
phy_wTree <- merge_phyloseq (phydat, tree)

## calculate distances for each method; group into single data.frame for plot
dist_bc <- phyloseq::distance(phy_wTree, "bray", binary = FALSE)
pcoa_bc <- ordinate(phy_wTree, method = "PCoA", distance = dist_bc)
dat_bc <- data.frame(pcoa_bc$vectors[,1:3]) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="bc")

dist_ds <- phyloseq::distance(phy_wTree, "bray", binary = TRUE)
pcoa_ds <- ordinate(phy_wTree, method = "PCoA", distance = dist_ds)
dat_ds <- data.frame(pcoa_ds$vectors[,1:3]) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="ds")

dist_mh <- phyloseq::distance(phy_wTree, "morisita", binary = FALSE)
pcoa_mh <- ordinate(phy_wTree, method = "PCoA", distance = dist_mh)
dat_mh <- data.frame(pcoa_mh$vectors[,1:3]) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="mh")

pcoa_wu = ordinate(phy_wTree, method="PCoA", distance="unifrac", weighted=TRUE)
dat_wu <- data.frame(pcoa_wu$vectors[,1:3]) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="wu")

pcoa_uu = ordinate(phy_wTree, method="PCoA", distance="unifrac", weighted=FALSE) ## note the Phyloseq impelentation normalizes branch length values
dat_uu <- data.frame(pcoa_uu$vectors[,1:3]) %>% mutate(SampleID = row.names(.)) %>% mutate(Measure="uu")

bigplot_df <- rbind(dat_bc, dat_ds, dat_mh, dat_wu, dat_uu)
rm(dat_bc, dat_ds, dat_mh, dat_wu, dat_uu)
bigplot_df <- merge(bigplot_df, meta)

## plot setup
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)
bigplot_df$Month <- factor(bigplot_df$Month, levels=c("June", "July", "September"))
bigplot_df$Measure <- factor(bigplot_df$Measure, levels=c('ds', 'bc', 'mh', 'uu', 'wu'))

## save plot as 'pcoa_fiveMetric_onlyGuano'; export at 1000x817
## note this plot isn't doing as good a job as the individual plots because:
## ... it's not showing the proportion of variation on each PC (which varies by plot)
## ... so maybe better not to facet except to use as a '30,000 foot' view to show broad trends at once
ggplot(bigplot_df, aes(x=Axis.1, y=Axis.2, color=Month, shape=Site)) +
  geom_point() +
  scale_color_manual(values=v3pal) +
  facet_wrap(Measure ~ ., nrow=2) +
  theme_devon() +
  theme(legend.position = "top")


## plotting individually directly with phyloseq is faster and better, with a bit of ggplot help
## going to add in the Egienvectors for Axis1/2 for each plot here
labaxisnames <- function(data, distance_type){
  axis1_name = data$values$Eigenvalues[1]
  axis2_name = data$values$Eigenvalues[2]
  tmp_df = data.frame(axis1_name, axis2_name)
  tmp_df %>% mutate(Measure=distance_type)
}

ds_name <- labaxisnames(pcoa_ds, "ds")
bc_name <- labaxisnames(pcoa_bc, "bc")
mh_name <- labaxisnames(pcoa_mh, "mh")
uu_name <- labaxisnames(pcoa_uu, "uu")
wu_name <- labaxisnames(pcoa_wu, "wu")

axis_names <- rbind(ds_name, bc_name, mh_name, uu_name, wu_name)

bigplot_df <- merge(bigplot_df, axis_names)

## gather each dataset
plot_bc <- plot_ordination(phy_wTree, pcoa_bc, color = "Month", shape = "Site")
plot_bc$data$Month <- factor(plot_bc$data$Month, levels = c("June", "July", "September"))

## plot; save as 'pcoa_bc_wellipse'; export at 650x550
pbc <- plot_bc + 
  geom_point(size = 4, alpha=0.8) + 
  scale_color_manual(values=v3pal) +
  theme_devon() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(subtitle="Bray-Curtis") +
  stat_ellipse(aes(group=Month))
#labs(subtitle="Bray-Curtis distance estimate") +
#stat_ellipse(data = plot_bc$data %>% filter(Site == "EN"), linetype = "dashed") + ## ellipse for EN months
#stat_ellipse(data = plot_bc$data %>% filter(Site == "HB"), linetype = "solid")  ## ellips for HB months
#stat_ellipse(aes(group = Month), linetype="longdash") ## this makes centroids just ~Month
#stat_ellipse(aes(group = Site), linetype="solid", color="black")
rm(plot_bc)

## plot; save as 'pcoa_ds_wellipse'; export at 650x550
plot_ds <- plot_ordination(phy_wTree, pcoa_ds, color = "Month", shape = "Site")
plot_ds$data$Month <- factor(plot_ds$data$Month, levels = c("June", "July", "September"))

pds <- plot_ds + 
  geom_point(size = 4, alpha=0.8) + 
  scale_color_manual(values=v3pal) +
  theme_devon() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(subtitle="Dice-Sorensen") +
  stat_ellipse(aes(group=Month))
rm(plot_ds)

## plot; save as 'pcoa_mh_wellipse'; export at 650x550
plot_mh <- plot_ordination(phy_wTree, pcoa_mh, color = "Month", shape = "Site")
plot_mh$data$Month <- factor(plot_mh$data$Month, levels = c("June", "July", "September"))
pmh <- plot_mh + 
  geom_point(size = 4, alpha=0.8) + 
  scale_color_manual(values=v3pal) +
  theme_devon() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(subtitle="Morisita-Horn") +
  stat_ellipse(aes(group=Month))
rm(plot_mh)


## plot; save as 'pcoa_uu_wellipse'; export at 650x550
plot_uu <- plot_ordination(phy_wTree, pcoa_uu, color = "Month", shape = "Site")
plot_uu$data$Month <- factor(plot_uu$data$Month, levels = c("June", "July", "September"))
puu <- plot_uu + 
  geom_point(size = 4, alpha=0.8) + 
  scale_color_manual(values=v3pal) +
  theme_devon() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(subtitle="Unweighted Unifrac") +
  stat_ellipse(aes(group=Month))
rm(plot_uu)

## plot; save as 'pcoa_wu_wellipse'; export at 650x550
plot_wu <- plot_ordination(phy_wTree, pcoa_wu, color = "Month", shape = "Site")
plot_wu$data$Month <- factor(plot_wu$data$Month, levels = c("June", "July", "September"))
pwu <- plot_wu + 
  geom_point(size = 4, alpha=0.8) + 
  scale_color_manual(values=v3pal) +
  theme_devon() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right") +
  labs(subtitle="Weighted Unifrac") +
  stat_ellipse(aes(group=Month))

rm(plot_wu)


### can plot all of these five with their unique axes % variance collectively:
require(cowplot)
top_row <- plot_grid(pds, pbc, pmh, nrow=1, labels = c("A", "B", "C"))
bottom_row <- plot_grid(puu, pwu, NULL, nrow = 1, rel_widths = c(1, 1.3, .7), labels = c("D", "E"))
## plot; save as 'pcoa_fiveMetric_wEigVal_wEllipse'; export at 1200x800
plot_grid(top_row, bottom_row, nrow=2)
rm(pds, pbc, pmh, puu, pwu)

########################################################################
## Part 3 == Adonis (multi factorial PERMANOVA) testing Site and Month main effects
##        == Testing if there are differences between centroids among groups
##        == first of two statistical tests; see PERMDISP in next section too
## for PERMANOVA details see: https://onlinelibrary.wiley.com/doi/full/10.1002/9781118445112.stat07841
########################################################################

## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method
dist_bc <- phyloseq::distance(phy_wTree, "bray")
dist_ds <- phyloseq::distance(phy_wTree, "bray", binary = TRUE)
dist_mh <- phyloseq::distance(phy_wTree, "morisita", binary = FALSE)
dist_uu <- phyloseq::distance(phy_wTree, "unifrac", weighted=FALSE)
dist_wu <- phyloseq::distance(phy_wTree, "wunifrac")

## export metadata from physeq object; do this because you've dropped a few samples for the ordination...
tmp.meta <- data.frame(sample_data(phy_wTree))

## run Adonis per distance method
adonis_ds <- adonis(dist_ds ~ Month * Site, data = tmp.meta) 
capture.output(data.frame(adonis_ds$aov.tab), file = "~/Repos/mysosoup/data/text_tables/adonis/adonis_ds.txt")

adonis_bc <- adonis(dist_bc ~ Month * Site, data = tmp.meta)
capture.output(data.frame(adonis_bc$aov.tab), file = "~/Repos/mysosoup/data/text_tables/adonis/adonis_bc.txt")

adonis_mh <- adonis(dist_mh ~ Month * Site, data = tmp.meta) 
capture.output(data.frame(adonis_mh$aov.tab), file = "~/Repos/mysosoup/data/text_tables/adonis/adonis_mh.txt")

adonis_uu <- adonis(dist_uu ~ Month * Site, data = tmp.meta) 
capture.output(data.frame(adonis_uu$aov.tab), file = "~/Repos/mysosoup/data/text_tables/adonis/adonis_uu.txt")

adonis_wu <- adonis(dist_wu ~ Month * Site, data = tmp.meta) 
capture.output(data.frame(adonis_wu$aov.tab), file = "~/Repos/mysosoup/data/text_tables/adonis/adonis_wu.txt")


########################################################################
## Part 4 == PERMDISP (multi factorial PERMANOVA) testing if within-group distances to group centroid differ across groups
##        == p < 0.xx result here complicates interpretation of Adoins (groups might be different because ..
##        == .. their within-group dispersions are variable among groups (not just differences among group centroids)
## For dispersions discussion, see: https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/12-2010.1 
########################################################################

## run permdisper (betadisper) for Site, then for Month groups for each distance measure (don't combine the term!):
ds_disper_Site <- betadisper(d = dist_ds, group =  tmp.meta$Site, type = c("median"))
ds_disper_Month <- betadisper(d = dist_ds, group =  tmp.meta$Month, type = c("median"))
capture.output(anova(ds_disper_Site), file="~/Repos/mysosoup/data/text_tables/permdisp/ds_Site_disper.txt")
capture.output(anova(ds_disper_Month), file="~/Repos/mysosoup/data/text_tables/permdisp/ds_Month_disper.txt")
##notrun: plot(ds_disper_Site)
##notrun: plot(ds_disper_Month)

bc_disper_Site <- betadisper(d = dist_bc, group =  tmp.meta$Site, type = c("median"))
bc_disper_Month <- betadisper(d = dist_bc, group =  tmp.meta$Month, type = c("median"))
capture.output(anova(bc_disper_Site), file="~/Repos/mysosoup/data/text_tables/permdisp/bc_Site_disper.txt")
capture.output(anova(bc_disper_Month), file="~/Repos/mysosoup/data/text_tables/permdisp/bc_Month_disper.txt")
##notrun: plot(bc_disper_Site)
##notrun: plot(bc_disper_Month)

mh_disper_Site <- betadisper(d = dist_mh, group =  tmp.meta$Site, type = c("median"))
mh_disper_Month <- betadisper(d = dist_mh, group =  tmp.meta$Month, type = c("median"))
capture.output(anova(mh_disper_Site), file="~/Repos/mysosoup/data/text_tables/permdisp/mh_Site_disper.txt")
capture.output(anova(mh_disper_Month), file="~/Repos/mysosoup/data/text_tables/permdisp/mh_Month_disper.txt")
##notrun: plot(mh_disper_Site)
##notrun: plot(mh_disper_Month)

uu_disper_Site <- betadisper(d = dist_uu, group =  tmp.meta$Site, type = c("median"))
uu_disper_Month <- betadisper(d = dist_uu, group =  tmp.meta$Month, type = c("median"))
capture.output(anova(uu_disper_Site), file="~/Repos/mysosoup/data/text_tables/permdisp/uu_Site_disper.txt")
capture.output(anova(uu_disper_Month), file="~/Repos/mysosoup/data/text_tables/permdisp/uu_Month_disper.txt")
##notrun: plot(uu_disper_Site)
##notrun: plot(uu_disper_Month)

wu_disper_Site <- betadisper(d = dist_wu, group =  tmp.meta$Site, type = c("median"))
wu_disper_Month <- betadisper(d = dist_wu, group =  tmp.meta$Month, type = c("median"))
capture.output(anova(wu_disper_Site), file="~/Repos/mysosoup/data/text_tables/permdisp/wu_Site_disper.txt")
capture.output(anova(wu_disper_Month), file="~/Repos/mysosoup/data/text_tables/permdisp/wu_Month_disper.txt")
##notrun: plot(wu_disper_Site)
##notrun: plot(wu_disper_Month)

## finally running Tukey's HSD on each dataset to identify any pairwise differences that were significant 
## doing so only for Month (given that there are only 2 levels for Site factor)
ds_month_hsd <- TukeyHSD(ds_disper_Month)
capture.output(data.frame(ds_month_hsd[[1]]) %>% mutate(Pairs=row.names(.)), file="~/Repos/mysosoup/data/text_tables/permdisp/ds_Month_tukeyHSD.txt")
bc_month_hsd <- TukeyHSD(bc_disper_Month)
capture.output(data.frame(bc_month_hsd[[1]]) %>% mutate(Pairs=row.names(.)), file="~/Repos/mysosoup/data/text_tables/permdisp/bc_Month_tukeyHSD.txt")
mh_month_hsd <- TukeyHSD(mh_disper_Month)
capture.output(data.frame(mh_month_hsd[[1]]) %>% mutate(Pairs=row.names(.)), file="~/Repos/mysosoup/data/text_tables/permdisp/mh_Month_tukeyHSD.txt")
uu_month_hsd <- TukeyHSD(uu_disper_Month)
capture.output(data.frame(uu_month_hsd[[1]]) %>% mutate(Pairs=row.names(.)), file="~/Repos/mysosoup/data/text_tables/permdisp/uu_Month_tukeyHSD.txt")
wu_month_hsd <- TukeyHSD(wu_disper_Month)
capture.output(data.frame(wu_month_hsd[[1]]) %>% mutate(Pairs=row.names(.)), file="~/Repos/mysosoup/data/text_tables/permdisp/wu_Month_tukeyHSD.txt")
