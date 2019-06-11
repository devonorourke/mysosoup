library(tidyverse)
library(qiime2R)
library(vegan)
library(viridis)


########################################################################
## section 1 - just the basics: generating anovas, plotting in base R
########################################################################

## import ASV table from qza:
tmpqza <- read_qza("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza")

## import metadata 
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv.gz", col_names = TRUE)
meta <- meta %>% 
  select(SampleID, Roost, CollectionMonth, SampleType, Site) %>% 
  rename('Month' = CollectionMonth) %>% 
  filter(SampleType == 'sample')
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$Month <- ifelse(meta$Month == "6", gsub("6", "June", meta$Month), meta$Month)
meta$Month <- ifelse(meta$Month == "7", gsub("7", "July", meta$Month), meta$Month)
meta$Month <- ifelse(meta$Month == "9", gsub("9", "September", meta$Month), meta$Month)
meta$SiteMonth <- paste(meta$Site, meta$Month, sep="-")

## retain only meta data sampleIDs remaining in table:
meta <- data.frame(meta) %>% filter(SampleID %in% colnames(tmpqza$data))

## run adonis for each data type
ds_adonis <- adonis(t(tmpqza$data) ~ Month * Site, data = meta, method = "bray", binary=TRUE)
ds_adonis_df <- as.data.frame(ds_adonis$aov.tab)
bc_adonis <- adonis(t(tmpqza$data) ~ Month * Site, data = meta, method = "bray", binary=FALSE)
bc_adonis_df <- as.data.frame(bc_adonis$aov.tab)

  ## because both main effects and interactions are significant, this indicates that either:
    ## 1. the centroid (or spatial median) among groups is different, or...
    ## 2. the within-group dispersion among groups is different
  ## we therefore need to test the second effect (within-group dispersions)

## run permdisper (betadisper) for Site, then for Month groups for each distance measure (don't combine the term!):
ds_disper_Site <- betadisper(d = vegdist(t(tmpqza$data), method = "bray", binary=TRUE), group =  meta$Site, type = c("median"))
ds_disper_Month <- betadisper(d = vegdist(t(tmpqza$data), method = "bray", binary=TRUE), group =  meta$Month, type = c("median"))
anova(ds_disper_Site) ## significant
anova(ds_disper_Month)  ## significant

bc_disper_Site <- betadisper(d = vegdist(t(tmpqza$data), method = "bray", binary=FALSE), group =  meta$Site, type = c("median"))
bc_disper_Month <- betadisper(d = vegdist(t(tmpqza$data), method = "bray", binary=FALSE), group =  meta$Month, type = c("median"))
anova(bc_disper_Site)   ## still significant, but a few orders of magnitude less
anova(bc_disper_Month)  ## still significant, but a few orders of magnitude less
  ## both distances are significant... so, let's plot to identify if adonis was mostly significant based on:
    ## 1. centroid
    ## 2. within-group dispersions

## quick plots; ellipse represent 1 sd around group medians; hull represents outermost points in 2d space
plot(ds_disper_Month, label = T, hull = T, ellipse = T)
plot(ds_disper_Site, label = T, hull = T, ellipse = T)
plot(bc_disper_Month, label = T, hull = T, ellipse = T)
plot(bc_disper_Site, label = T, hull = T, ellipse = T)

    ## Summary: it looks like occurrence data partitions differences between groups much more than abundance info ...
    ## This may likely be because abundance data is swamping out low abundant differences (it's not clear though if low abundance differences are what differentiate the groups!)


boxplot(ds_disper_Site) ## demonstrates how within-group dispersions vary (what permdisp/betadisp is testing!!)
boxplot(ds_disper_Month)
boxplot(bc_disper_Site) ## more dispersion at EN site
boxplot(bc_disper_Month)  ## more dispersion in June

## which pairwise differences are most significant? (not testing for Site, given that there is only one pair!)
ds_hsd_month <- TukeyHSD(ds_disper_Month)
ds_hsd_month_df <- data.frame(ds_hsd_month[[1]]) %>% mutate(Pairs=row.names(.))   ## June:July and June:September dispersions differ

bc_hsd_month <- TukeyHSD(bc_disper_Month)
bc_hsd_month_df <- data.frame(bc_hsd_month[[1]]) %>% mutate(Pairs=row.names(.))   ## just June:September differences here


########################################################################
## section 2 - enhanced plots - not that useful...
########################################################################

## function for plot theme:
theme_devon <- function () { theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(panel.background  = element_blank(), plot.background = element_rect(fill="transparent", colour=NA), 
          legend.background = element_rect(fill="transparent", colour=NA), legend.key = element_rect(fill="transparent", colour=NA))}

## function to make custom dispersion plots
dispplotdat.function <- function(displist, DistMethod, GroupVar){
  tmp.disp_sco <- as.data.frame(scores(displist, display = "sites"))
  tmp.disp_sco$bigID <- row.names(tmp.disp_sco)
  row.names(tmp.disp_sco) <- NULL
  colnames(tmp.disp_sco) <- c("pc1end", "pc2end", "SampleID")
  tmp.disp_sco <- merge(tmp.disp_sco, meta) %>% select(-Roost, -SampleType)
  tmp.disp_centroids <- as.data.frame(scores(displist, display="centroids"))  
  tmp.disp_centroids$SiteMonth <- row.names(tmp.disp_centroids)
  row.names(tmp.disp_centroids) <- NULL
  colnames(tmp.disp_centroids) <- c("pc1center", "pc2center", GroupVar)
  tmp.disp.plot <- merge(tmp.disp_sco, tmp.disp_centroids)  
  tmp.disp.plot <- tmp.disp.plot %>% mutate(DistMethod=DistMethod)
}


## generate point coordinates:
ds_site_plot_df <- dispplotdat.function(ds_disper_Site, "ds", "Site")

## generate fraction Eigenval (proportion variation)
tmplabs <- data.frame(round(100*ds_disper_Site$eig / sum(ds_disper_Site$eig), 2)) %>% mutate(PCaxis = row.names(.))
colnames(tmplabs)[1] <- 'EigVal'
tmpxval <- tmplabs %>% filter(PCaxis == 'PCoA1') %>% select(EigVal) %>% pull()
tmpyval <- tmplabs %>% filter(PCaxis == 'PCoA2') %>% select(EigVal) %>% pull()

## set levels for both objects:
ds_site_plot_df$Site <- factor(ds_site_plot_df$Site, levels=c("EN", "HB"))

## set palette:
vpal2 <- c("dodgerblue4", "firebrick")
## unused: vpal6 <- c("dodgerblue4", "firebrick", "darkorchid", "darkorange2", "darkslategray3", "goldenrod")

## plot
ggplot() +
  scale_color_manual(values=vpal2) +
  geom_segment(data=ds_site_plot_df, aes(x=pc1end, y=pc2end, xend=pc1center, yend=pc2center, color=Site), alpha=0.1) +
  geom_point(data=ds_site_plot_df, aes(x=pc1end, y=pc2end, color=Site), size=2, alpha=0.8) +
  theme_devon() +
  theme(legend.position = 'top') +
  guides(color=guide_legend(nrow = 2)) +
  labs(x=paste("PCoA-1 (", tmpxval, " %)"),
       y=paste("PCoA-2 (", tmpyval, " %)"),
       color="", fill="")


### other unused code:
## generate hull coordinates
#   library(plyr)
#   find_hull <- function(df) df[chull(df$pc1end, df$pc2end), ]
#   hulls <- ddply(tmp, "SiteMonth", find_hull)


# plot
#   ggplot() + scale_color_manual(values=vpal6) +
#  scale_fill_manual(values=vpal6) +
#  geom_polygon(data=hulls, aes(x=pc1end, y=pc2end, color=SiteMonth), fill="white", alpha=0.01) +
#  geom_segment(data=tmp, aes(x=pc1end, y=pc2end, xend=pc1center, yend=pc2center, color=SiteMonth), alpha=0.1) +
#  geom_point(data=tmp, aes(x=pc1end, y=pc2end, color=SiteMonth), size=2, alpha=0.8) +
#  stat_ellipse(data = tmp, aes(x=pc1end, y=pc2end, color=SiteMonth), linetype = "dashed", size=0.5, type="t") +
#  theme_devon() +
#  theme(legend.position = 'top') +
#  guides(color=guide_legend(nrow = 2)) + labs(x=paste("PCoA-1 (", tmpxval, " %)"), y=paste("PCoA-2 (", tmpyval, " %)"), color="", fill="")