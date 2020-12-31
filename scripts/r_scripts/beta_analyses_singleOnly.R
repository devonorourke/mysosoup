library(tidyverse)
library(qiime2R)
library(reshape2)
library(phyloseq)
library(ape)
library(vegan)
library(ggpubr)
library(svglite)

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
## Part 0 - data import: metadata, sequence abundance .qza object, and newick tree
########################################################################

## import metadata
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE)


## add taxonomy information
## amend import of taxa info
taxa <- read_csv(file="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/taxonomy/filtd_tax_dataframe_ALL.csv")
colnames(taxa)[1] <- "ASVid"

## import list of nonMYSO samples: 
nonMYSO_sampleID_list <- read_table(file="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/host/nonMYSO_sampleIDlist.txt", col_names = FALSE)

## import OTU table from QZA artifact
## download qza file:
# download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/Mangan.clust_p985_table_Rarefyd.qza", "tmp.qza")
# then set PATH to whatever directory needed:
# qzapath = "PATH/TO/tmp.qza"

## import read abundance information from local .qza file...
## set your own path below!
## filter out pooled samples and analyze only single samples
qzapath = "/Users/devonorourke/github/mysosoup/data/qiime_qza/Mangan.clust_p985_table_Rarefyd.qza"  ## amend this line to your own path!
featuretable <- read_qza(qzapath)
mat.tmp <- featuretable$data
rm(featuretable)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
rm(df.tmp)
colnames(tmp) <- c("ASVid", "SampleID", "Reads")  ## note these are OTUs, but the original taxonomy was assigned to individual ASVs
df.tmp <- merge(tmp, taxa)
rm(tmp)
df.tmp <- merge(df.tmp, meta)
df.tmp <- df.tmp %>% filter(BatchType == "single" & SampleType == "sample" & !SampleID %in% nonMYSO_sampleID_list$X1)
mat.tmp <- dcast(data = df.tmp, formula = SampleID ~ ASVid, value.var='Reads', fill = 0)
row.names(mat.tmp) <- mat.tmp$SampleID
mat.tmp$SampleID <- NULL
rm(df.tmp, qzapath)

## create phyloseq object with necessary samples/taxa/meta info
OTU <- otu_table(mat.tmp, taxa_are_rows = FALSE) ## import as physeq object
meta <- meta %>% filter(SampleID %in% row.names(mat.tmp))
row.names(meta) <- meta$SampleID
row.names(taxa) <- taxa$ASVid
phy_tax <- as.matrix(taxa) %>% tax_table(.)
phydat <- phyloseq(OTU, phy_tax, meta)

rm(mat.tmp, phy_tax, OTU)

## sanity check:
nsamples(phydat)  ## 189 expected
ntaxa(phydat) ## 1057 taxa
sample_names(phydat)  ## sample names as SampleID in meta
  ## all good.


## import tree object to for phylogenetic distance measures
## download qza file:
# download.file("https://github.com/devonorourke/mysosoup/raw/master/data/trees/Mangan.clust_p985_rooted_tree.nwk", "tree.nwk")
# tree <- read.tree(file = "~/Desktop/tree.nwk")

## alternatively run from local:
tree <- read.tree(file = "~/github/mysosoup/data/trees/Mangan.clust_p985_rooted_tree.nwk")

## add tree info to physeq object
phy_wTree <- merge_phyloseq (phydat, tree)


########################################################################
## Part 1a: calculate distances using one of meach metric:
## abundance unweighted, phylogenetic unweighted (Dice-Sorenson)
## abundance weighted, phylogenetic unweighted (Bray-Curtis)
## abundance unweighted, phylogenetic weighted (Unifrac unweighted)
## abundance weighted, phylogenetic weighted (Unifrac weighted)

## Part 1b: make the plots of each PCoA
########################################################################

## calculate distances for each method; group into single data.frame for plot
beta_nonphylo_plotdat_function <- function(distance_method, distance_alias, binary_choice){
  dist_tmp <- phyloseq::distance(phy_wTree, distance_method, binary = binary_choice)
  pcoa_tmp <- ordinate(phy_wTree, method = "PCoA", distance = dist_tmp)
  Props = pcoa_tmp$values$Eigenvalues/sum(pcoa_tmp$values$Eigenvalues)
  Prop1 = Props[1]
  Prop2 = Props[2]
  data.frame(pcoa_tmp$vectors[,1:2]) %>% 
    mutate(SampleID = row.names(.)) %>% 
    mutate(Measure=distance_alias) %>% 
    mutate(Axis1var = Prop1, Axis2var = Prop2)
}

dat_ds <- beta_nonphylo_plotdat_function("bray", "Dice-Sorensen", TRUE)
dat_bc <- beta_nonphylo_plotdat_function("bray", "Bray-Curtis", FALSE)


beta_phylo_plotdat_function <- function(distance_alias, weight_choice){
  pcoa_tmp <- ordinate(phy_wTree, method = "PCoA", distance = "unifrac", weighted = weight_choice)
  Props = pcoa_tmp$values$Eigenvalues/sum(pcoa_tmp$values$Eigenvalues)
  Prop1 = Props[1]
  Prop2 = Props[2]
  data.frame(pcoa_tmp$vectors[,1:2]) %>% 
    mutate(SampleID = row.names(.)) %>% 
    mutate(Measure=distance_alias) %>% 
    mutate(Axis1var = Prop1, Axis2var = Prop2)
}

dat_uu <- beta_phylo_plotdat_function("Unifrac Unweighted", FALSE)
dat_wu <- beta_phylo_plotdat_function("Unifrac Weighted", TRUE)

## combine the dataset for a single faceted plot
bigplot_df <- rbind(dat_bc, dat_ds, dat_wu, dat_uu)
bigplot_df <- merge(bigplot_df, meta)
rm(dat_bc, dat_ds, dat_wu, dat_uu)

bigplot_df <- bigplot_df %>% 
  mutate(MonthVal = case_when(CollectionMonth == "6" ~ "June",
                              CollectionMonth == "7" ~ "July",
                              CollectionMonth == "9" ~ "Sept"),
         SiteVal = case_when(Site == "Egner" ~ "EN",
                             Site == "HickoryBottoms" ~ "HB"))

##### 1b. the plot:
## plot setup
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)
bigplot_df$MonthVal <- factor(bigplot_df$MonthVal, levels=c("June", "July", "Sept"))
bigplot_df$Measure <- factor(bigplot_df$Measure, levels=c('Dice-Sorensen', 'Bray-Curtis', 'Unifrac Unweighted', 'Unifrac Weighted'))
bigplot_df$Axis1var <- as.character(bigplot_df$Axis1var)
bigplot_df$Axis2var <- as.character(bigplot_df$Axis2var)

## max variance for PCs 1 and 2 across all plots (to keep coordinate distances similar)
data.frame(maxX = max(bigplot_df$Axis.1),
           minX = min(bigplot_df$Axis.1),
           maxY = max(bigplot_df$Axis.2),
           minY = min(bigplot_df$Axis.2))
## use x = c(-.5, 0.4)
## use y = c(-.4, .5)

## create function to individually plot each % variance explained by PC1/2 inset on bottom left of each facet
pcoaPlotfunction <- function(BetaMetric){
  pc1textval = bigplot_df %>% filter(Measure == BetaMetric) %>% distinct(Axis1var) %>% pull() %>% as.numeric(.) %>% round(3)*100
  pc2textval = bigplot_df %>% filter(Measure == BetaMetric) %>% distinct(Axis2var) %>% pull() %>% as.numeric(.) %>% round(3)*100
  ggplot(bigplot_df %>% filter(Measure == BetaMetric), 
         aes(x=Axis.1, y=Axis.2, color=MonthVal, shape=SiteVal)) +
    geom_point() +
    scale_color_manual(values=v3pal) +
    facet_wrap(Measure ~ ., nrow=2) +
    theme_devon() +
    theme(strip.text = element_text(size=11)) +
    lims(x=c(-0.5, 0.4),
         y=c(-0.4, 0.5)) +
    labs(color="Month", shape="Site", 
         x=paste0("\nPC1:  ", pc1textval, "%\n"), 
         y=paste0("PC2:  ", pc2textval, "%\n"))
}


## make individual plots
## save individual plots in case of need to add individually in Illustrator later
p_ds <- pcoaPlotfunction("Dice-Sorensen")
p_bc <- pcoaPlotfunction("Bray-Curtis")
p_uu <- pcoaPlotfunction("Unifrac Unweighted")
p_wu <- pcoaPlotfunction("Unifrac Weighted")

p_ds
ggsave("~/github/mysosoup/figures/figure_3s_ds_pcoa.png", width = 15, height = 15, units = "cm")
ggsave("~/github/mysosoup/figures/figure_3s_ds_pcoa.svg", width = 15, height = 15, units = "cm")

p_bc
ggsave("~/github/mysosoup/figures/figure_3s_bc_pcoa.png", width = 15, height = 15, units = "cm")
ggsave("~/github/mysosoup/figures/figure_3s_bc_pcoa.svg", width = 15, height = 15, units = "cm")

p_uu
ggsave("~/github/mysosoup/figures/figure_3s_uu_pcoa.png", width = 15, height = 15, units = "cm")
ggsave("~/github/mysosoup/figures/figure_3s_uu_pcoa.svg", width = 15, height = 15, units = "cm")

p_wu
ggsave("~/github/mysosoup/figures/figure_3s_wu_pcoa.png", width = 15, height = 15, units = "cm")
ggsave("~/github/mysosoup/figures/figure_3s_wu_pcoa.svg", width = 15, height = 15, units = "cm")

## stitch together
ggarrange(p_ds, p_bc, p_uu, p_wu, nrow = 2, ncol = 2, common.legend = TRUE)
ggsave("~/github/mysosoup/figures/figure3_betaOrd.png", width = 20, height = 20, units = "cm")
ggsave("~/github/mysosoup/figures/figure3_betaOrd.svg", width = 20, height = 20, units = "cm")

rm(p_ds, p_bc, p_uu, p_wu)

########################################################################
## Part 2 == Adonis (multi factorial PERMANOVA) testing Site and Month main effects
##        == Testing if there are differences between centroids among groups
##        == first of two statistical tests; see PERMDISP in next section too
## for PERMANOVA details see: https://onlinelibrary.wiley.com/doi/full/10.1002/9781118445112.stat07841
########################################################################

## calculate distances individually, then run PERMANOVA (via Adonis in Vegan); one per distance method
dist_ds <- phyloseq::distance(phy_wTree, "bray", binary = TRUE)
dist_bc <- phyloseq::distance(phy_wTree, "bray")
dist_uu <- phyloseq::distance(phy_wTree, "unifrac", weighted=FALSE)
dist_wu <- phyloseq::distance(phy_wTree, "wunifrac")

## run Adonis per distance method
adonis_data_function <- function(distanceData, distanceMetric){
  adonis_tmp <- adonis(distanceData ~ CollectionMonth * Site, data = meta)
  data.frame(adonis_tmp$aov.tab) %>% mutate(Metric = distanceMetric,
                                            Class = row.names(.))
}

adonis_ds <- adonis_data_function(dist_ds, "Dice-Sorensen")
adonis_bc <- adonis_data_function(dist_bc, "Bray-Curtis")
adonis_uu <- adonis_data_function(dist_uu, "Unifrac Unweighted")
adonis_wu <- adonis_data_function(dist_uu, "Unifrac Weighted")

adonis_all <- rbind(adonis_ds, adonis_bc, adonis_uu, adonis_wu)
adonis_all <- adonis_all[,c(8,1,2,3,4,5,6,7)]
write.csv(adonis_all, quote=FALSE, row.names = FALSE,
          file = "~/github/mysosoup/data/text_tables/adonis/adonis_data_allMetrics.csv")

rm(adonis_ds, adonis_bc, adonis_uu, adonis_wu)
########################################################################
## Part 3 == PERMDISP (multi factorial PERMANOVA) testing if within-group distances to group centroid differ across groups
##        == p < 0.xx result here complicates interpretation of Adoins (groups might be different because ..
##        == .. their within-group dispersions are variable among groups (not just differences among group centroids)
## For dispersions discussion, see: https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/12-2010.1 
########################################################################

## run permdisper (betadisper) for Site, then for Month groups for each distance measure (don't combine the term!):
beta_disper_function <- function(distanceData, distanceMetric){
  tmp_disper_Site_list <- betadisper(d = distanceData, group =  meta$Site, type = c("median"))
  tmp_disper_Site_anova <- data.frame(anova(tmp_disper_Site_list))
  tmp_disper_Site_anova <- tmp_disper_Site_anova %>% 
    mutate(Class = row.names(.)) %>% mutate(Metric = distanceMetric) %>% mutate(Factor = "Site")
  tmp_disper_Month_list <- betadisper(d = distanceData, group =  meta$CollectionMonth, type = c("median"))
  tmp_disper_Month_anova <- data.frame(anova(tmp_disper_Month_list))
  tmp_disper_Month_anova <- tmp_disper_Month_anova %>% 
    mutate(Class = row.names(.)) %>% mutate(Metric = distanceMetric) %>% mutate(Factor = "Month")
  tmp_out <- rbind(tmp_disper_Site_anova, tmp_disper_Month_anova)
  tmp_out[,c(6,1,2,3,4,5,8,7)]
}

bidsp_ds <- beta_disper_function(dist_ds, "Dice-Sorensen")
bidsp_bc <- beta_disper_function(dist_bc, "Bray-Curtis")
bidsp_uu <- beta_disper_function(dist_uu, "Unifrac Unweighted")
bidsp_wu <- beta_disper_function(dist_wu, "Unifrac Weighted")

bdisp_all <- rbind(bidsp_ds, bidsp_bc, bidsp_uu, bidsp_wu) %>% 
  mutate(Sum.Sq = round(Sum.Sq,4),
         Mean.Sq = round(Sum.Sq,4),
         F.value = round(Sum.Sq,4),
         Pr..F. = round(Pr..F., 4))

write.csv(bdisp_all, quote=FALSE, row.names = FALSE,
          file = "~/github/mysosoup/data/text_tables/permdisp/betaDispersion_data_allMetrics.csv")
