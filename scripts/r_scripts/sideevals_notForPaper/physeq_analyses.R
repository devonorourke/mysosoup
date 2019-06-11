### WARNING --- IT'S CLEAR THAT SOMETHING FUNNY IS GOING ON WITH A FEW SAMPLES
### TO DROP THESE SAMPLES IN PHYLOSEQ:
'%!in%' <- function(x,y)!('%in%'(x,y))
badsamples <- c('7272017EGC1', '7272017HBD3')
mrd = subset_samples(mr, SampleID %!in% badsamples)


library(tidyverse)
library(reshape2)
library(phyloseq)
library(vegan)
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


## import sequence data:
data <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/mangan.asvtable.long.csv.gz", col_names = TRUE)
data_mat <- dcast(data, SampleID ~ ASVid, value.var = 'Reads', fill = 0)
row.names(data_mat) <- data_mat$SampleID
NULL -> data_mat$SampleID
phy_otu <- otu_table(data_mat, taxa_are_rows = FALSE) ## import as physeq object 
#rm(data_mat, data)

## import metadata 
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/mangan_metadata.csv.gz", col_names = TRUE)
#meta <- meta %>% filter(SampleID %in% )
row.names(meta) <- meta$SampleID
phy_meta <- sample_data(meta)

## minor metadata table for some plots:
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
row.names(taxa) <- taxa$ASVid
NULL -> taxa$ASVid
tax_mat <- as.matrix(taxa)
phy_tax <- tax_table(tax_mat)
rm(tax_mat, taxa)


## combine as physeq object:
m <- phyloseq(phy_otu, phy_tax, phy_meta)
rm(phy_otu, phy_tax, phy_meta)

## rarefy:
mr = rarefy_even_depth(m, sample.size=9000, replace = FALSE, rngseed = 1000)
#smr = rarefy_even_depth(sm, sample.size=9000, replace = FALSE, rngseed = 1000)

## estimate richness
alpha_df <- estimate_richness(mr, split = TRUE, measures = c("Observed", "Simpson", "Shannon"))
alpha_df <- estimate_richness(smr, split = TRUE, measures = c("Observed", "Simpson", "Shannon"))
## plot alpha diversity
alpha_df$SampleID <- row.names(alpha_df)
alt_alpha <- alpha_df %>% 
  gather(data = ., "Observed", "Shannon", "Simpson", key = "Measure", value="value") %>%
  arrange(Measure, value)
alt_alpha$SampleID <- gsub("^X", "", alt_alpha$SampleID)
alt_alpha <- merge(alt_alpha, tinymeta)

ggplot(alt_alpha, aes(x=Labeler, y=value, color=CollectionMonth)) + 
  geom_jitter(width = 0.1, alpha=0.8) + 
  scale_color_manual(values=c("green", "orange", "red", "gray30")) +
  facet_wrap( ~ Measure, scales = 'free_y') +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1))


## estimate distances
vegdist.bray <- vegdist(data_mat, method = "bray", binary = FALSE)
vegdist.dice <- vegdist(data_mat, method = "bray", binary = TRUE)

## gather NMDS and plot
nmds.bray <- metaMDS(vegdist.bray, distance = "bray", k=3, autotransform = FALSE)
nmds.dice <- metaMDS(vegdist.dice, distance = "bray", k=3, autotransform = FALSE)

ggplot(data.frame(nmds.bray$points), aes(x=MDS1, y=MDS2)) + geom_point()
nmds.bray.df <- data.frame(nmds.bray$points)
nmds.bray.df$SampleID <- row.names(nmds.bray.df)
nmds.bray.df <- merge(nmds.bray.df, tinymeta)

## three axes: plotting 1 and 2
ggplot(nmds.bray.df, aes(x=MDS1, y=MDS2, color=CollectionMonth, shape=Site)) +
  geom_point(size = 2) +
  theme_devon()

## plotting 2 and 3
bray.mds23 <- ggplot(nmds.bray.df, aes(x=MDS2, y=MDS3, color=CollectionMonth, shape=Site)) +
  geom_point(size = 2) +
  theme_devon()

## generate tree:
m_tree = rtree(ntaxa(mr), rooted=TRUE, tip.label=taxa_names(mr))
plot_tree(m_tree, color="class_name", ladderize = "left") + coord_polar(theta = "y")
class(m_tree)
out <- merge_phyloseq(m_tree, mr)
plot_tree(out, color="class_name", ladderize = "left")

## run 
# ----------------------------------------------------------------


## standardize sampling depths to median value of total abundance per sample
total = median(sample_sums(m))
standf = function(x, t=total) round(t * (x / sum(x)))
ms = transform_sample_counts(m, standf)

## remove ASVs where total abundance of that ASV is less than 0.05% of mean ASV abundance (from normalized dataset):
msr  = transform_sample_counts(ms, function(x) x / sum(x) )
msrf = filter_taxa(msr, function(x) mean(x) > 1e-5, TRUE)
ntaxa(msrf)   ## 2620 ASVs remain (thus about half of all ASVs removed)
nsamples(msrf)  ## 297 samples
keepTaxa <- taxa_names(msrf)
mf = prune_taxa(keepTaxa, ms)  ## object is normalized with pruned ASVs (just 2620 remain)
ntaxa(mf) ## 2620 - samples remain

## examine alpha diversities
alpha_df <- estimate_richness(m)
alpha_df 
plot_richness(m)
## maybe the OTU tableis less than 0? look back at filter...



## examine beta diversities
## calculate distances:
raupDist <- distance(mf, method="raup")
brayDist <- distance(mf, method="bray")
moriDist <- distance(mf, method="morisita")

## which ordination method to choose: http://ordination.okstate.edu/overview.htm ?
# Calculate ordination via Vegan and NMDS
raupMDS  <- metaMDS(comm = raupDist, distance = "raup", k = 2, trymax = 30)
brayMDS  <- metaMDS(comm = brayDist, distance = "bray", k = 2, trymax = 30)
moriMDS  <- metaMDS(comm = moriDist, distance = "morisita", k = 2, trymax = 30)

# Combine NMDS results into single object
raup.df <- data.frame(raupMDS$points) %>% mutate(Dist="raup") %>% mutate(SampleID=row.names(raupMDS$points))
bray.df <- data.frame(brayMDS$points) %>% mutate(Dist="bray") %>% mutate(SampleID=row.names(brayMDS$points))
mori.df <- data.frame(moriMDS$points) %>% mutate(Dist="mori") %>% mutate(SampleID=row.names(moriMDS$points))
nmds.df <- rbind(raup.df, bray.df, mori.df)
rm(raup.df, bray.df, mori.df)

# add metadata:
nmds.df <- merge(nmds.df, meta)

## plot
ggplot(nmds.df, aes(x=MDS1, y=MDS2, color=ContamArea, label=SampleID)) + 
  geom_point() +
  geom_text(data=nmds.df %>% filter(MDS1 > 0.5)) +
  facet_grid(~ Dist)


# Calculate ordination via Phyloseq and DMS
p.raupMDS  <- ordinate(mf, method = "DCA", distance="raup")
p.brayMDS  <- ordinate(mf, method = "DCA", distance="bray")
p.moriMDS  <- ordinate(mf, method = "DCA", distance="morisita")

p.raup_df <- data.frame(p.raupMDS$rproj) %>% mutate(SampleID=row.names(.)) %>% mutate(Dist="raup")
p.bray_df <- data.frame(p.brayMDS$rproj) %>% mutate(SampleID=row.names(.)) %>% mutate(Dist="bray")
p.mori_df <- data.frame(p.moriMDS$rproj) %>% mutate(SampleID=row.names(.)) %>% mutate(Dist="mori")
p.nmds.df <- rbind(p.raup_df, p.bray_df, p.mori_df)
rm(p.raup_df, p.bray_df, p.mori_df)
p.nmds.df <- merge(p.nmds.df, meta)

## plot
ggplot(p.nmds.df, aes(x=DCA1, y=DCA2)) + geom_point() + facet_grid(~Dist)

##  --- phyloseq accessor functions --- ##
ntaxa(mang_phy) ## 5290 ASVs
nsamples(mang_phy)  ## 297 samples
sample_names(mang_phy)[1:5] ## first 5 sample names
rank_names(mang_phy)  ## possible taxa rank names to use in plots
sample_variables(mang_phy)  ## possible variables to filter/use in plots/calculations
otu_table(mang_phy)[1:5, 1:5] ## access {dim} of OTU table
tax_table(mang_phy)[1:5, 1:4] ## access {din} of TAX table


rank_names(m)

