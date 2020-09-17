library(tidyverse)
library(vegan)
library(scales)
library(qiime2R)
library(reshape2)
library(formattable)
library(ggpubr)
library(Matrix)
library(multcompView)
library(scico)
library(FSA)

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

################################################################################
## Step 0 - import metadata, taxonomy infomration, and abundance tables
################################################################################

## import metadata 
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE)
meta <- meta %>% select(SampleID, Roost, CollectionMonth, Site, SampleType, BatchType)
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$CollectionMonth[is.na(meta$CollectionMonth)] <- "control"
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "6", gsub("6", "June", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "7", gsub("7", "July", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "9", gsub("9", "September", meta$CollectionMonth), meta$CollectionMonth)
meta$Labeler <- paste(meta$Site, meta$Roost, sep="-")
meta$Labeler <- ifelse(meta$Labeler == "control-control", gsub("control-control", "control", meta$Labeler), meta$Labeler)
meta$CollectionMonth <- as.factor(meta$CollectionMonth)


## add taxonomy information
## amend import of taxa info
taxa <- read_csv(file="/scratch/dro49/qiimetmp/mysotmp/filtd_tax_dataframe_ALL.csv")
# taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
# taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
# taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
# taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
# taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
# taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"



## ..filter out remaining ASVs present exclusively in control samples, then calculate Hill Number values
## amend path to github file; need instruction on how to download first, then import from local
## run from local: qzapath = "~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza"
qzapath = "/scratch/dro49/qiimetmp/mysotmp/Mangan.clust_p985_table_Rarefyd.qza"
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
tmp.hill <- data.frame(renyi(mat.tmp, scales = c(0,1,2), hill=TRUE)) %>% mutate(SampleID = row.names(.))
colnames(tmp.hill)[1:3] <- c("Observed", "Shannons", "Simpsons")
Alpha_df <- gather(tmp.hill, key="Alpha_metric", value = "Alpha_value", c('Observed', 'Shannons', 'Simpsons'))
mangan_hill_df <- merge(Alpha_df, meta)
rm(tmp.hill, mat.tmp, Alpha_df, df.tmp, qzapath)


################################################################################
## Step 1 - run ANOVA for Site, CollectionMonth, and BatchType ...
## ... to identify main effects and interactions

## Analysis shows that BatchType is itself a main effect
## See associated boxplot 
################################################################################

#### Run ANOVA for group significance:
alpha_anova_allSamps <- function(qval, filename) {
  aov.out <- aov(Alpha_value ~ CollectionMonth * Site * BatchType, 
                 data=mangan_hill_df %>% filter(Alpha_metric == qval))
}

anova_observed_allSamps <- alpha_anova_allSamps("Observed")
summary(anova_observed_allSamps)
  ## sig Main effects: CollectionMonth, BatchType
  ## sig Interactions: Site:BatchType

anova_shannons_allSamps <- alpha_anova_allSamps("Shannons")
summary(anova_shannons_allSamps)
  ## sig Main effects: CollectionMonth, Site, BatchType
  ## sig Interactions: CollectionMonth:Site, Site:BatchType

anova_simpsons_allSamps <- alpha_anova_allSamps("Simpsons")
summary(anova_simpsons_allSamps)
  ## same as Shannons ^^

### plot to figure out effects of richness on 3 main effects
mangan_hill_df$CollectionMonth <- factor(mangan_hill_df$CollectionMonth, levels = c("June", "July", "September"))
#v3pal <- viridis::plasma(3, begin = 0.2, end = 0.95, direction = -1)
monthpal <- scico(3, palette = 'batlow', begin = 0.2, end = 0.8)

ggplot(mangan_hill_df, aes(x=Site, y=Alpha_value, color=CollectionMonth)) + 
  geom_boxplot(outlier.shape = NA, color="gray30") +
  #geom_jitter(width = 0.2, alpha=0.6) + 
  geom_point(position = position_jitterdodge()) +
  scale_color_manual(values = c(monthpal, "gray40"), labels=c("June", "July", "September")) +
  facet_grid(Alpha_metric ~ BatchType) +
  labs(x="\nCollection site", y="Estimated diversity", color = "Month") +
  theme_devon() + theme(legend.position = "top")

## Two effects are apparent between BatchTypes:
## 1. Hickory pooled samples have consistently highest diversity among all BatchType*Sites groups
## 2. Pooled samples have higher observed richness in both Sites

## Another point to notice is that we don't have sufficient sample sizes for pools if ...
## ... they are split by Site-Month. We also see a more bimodal distribution of diversity values ...
## ... among pooled samples, which given we didn't control the actual mass of guano might contribute to this variation

## Because BatchType is skewing these data and not a relevant main effect we're studying, ...
## ... we'll split up these analyses into 'single' and 'pooled' (we have sufficient samples sizes for both)

## cleanup
rm(anova_observed, anova_shannons, anova_simpsons, anovafunction)

################################################################################
## Step 2 - Rerun ANOVA for Site, CollectionMonth for EACH BatchType separately

## For SINGLE sample BatchTypes:
  ## we find significant effect for CollectionMonth and Site for Observed OTUs (richness), but...
  ## we find NO significant effects for Month or Site for Shannons/Simpsons

## For POOLED sample BatchTypes we observe different main effects for each alpha metric:
  ## significant effects for Site, but only weak effect of month (p=0.095)
  ## significant effects for Site and Month for Shannons and Simpsons (no interaction effects)

################################################################################

## First, examine Single BatchType samples
## update this function to export text table of summary outputs rather than printing to screen
alpha_anova_single <- function(qval, filename) {
  aov.out <- aov(Alpha_value ~ CollectionMonth * Site, 
                 data=mangan_hill_df %>% 
                   filter(BatchType == 'single') %>% 
                   filter(Alpha_metric == qval))
}

alpha_single_observed <- alpha_anova_single("Observed")
summary(alpha_single_observed)
  ## sig Main effects: CollectionMonth, Site
  ## no sig Interactions

alpha_single__shannons <- alpha_anova_single("Shannons")
summary(alpha_single__shannons)
  ## no sig Main effects; weak effect for CollectionMonth (p = 0.0713)
  ## no sig Interactions

alpha_single__simpsons <- alpha_anova_single("Simpsons")
summary(alpha_single__simpsons)
  ## no sig Main effects
  ## no sig Interactions


## Second, examine Pooled BatchType samples
alpha_anova_pool <- function(qval, filename) {
  aov.out <- aov(Alpha_value ~ CollectionMonth * Site, 
                 data=mangan_hill_df %>% 
                   filter(BatchType == 'pool') %>% 
                   filter(Alpha_metric == qval))
}

alpha_pool_observed <- alpha_anova_pool("Observed")
summary(alpha_pool_observed)
## sig Main effects: Site (weak effect for Month (p=0.095))
## no sig Interactions

alpha_pool_shannons <- alpha_anova_pool("Shannons")
summary(alpha_pool_shannons)
## sig Main effects for Month and Site
## no sig Interactions

alpha_pool_simpsons <- alpha_anova_pool("Simpsons")
summary(alpha_pool_simpsons)
## same as Shannons ^^

## these results suggest we should be analyzing the BatchType as separate experiments
## the question remains as to whether we should analyze these using a parametric test or not
## can also try running the test using Kruskal-Wallis?


################################################################################
## non parametric tests for ranked mean differences
## splitting samples of different BatchTypes (single vs. pooled samples) 

## Kruskal-Wallis tests are significant for all but single Shannon's diversity comparisons (for Site-Month groups)
## this matches our ANOVA analyses generally
################################################################################

###### Run Kruskal Wallis for nonparametric test of same data:
kwfunction <- function(qval, batch){
  mangan_hill_df$Grouper <- paste(mangan_hill_df$CollectionMonth, mangan_hill_df$Site, sep="-")
  mangan_hill_df$Grouper <- as.factor(mangan_hill_df$Grouper)
  kruskal.test(Alpha_value ~ Grouper, 
               data = mangan_hill_df %>% 
                 filter(Alpha_metric == qval) %>% 
                 filter(BatchType == batch))
}

kw_pool_observed <- kwfunction('Observed', 'pool')
  ## Kruskal-Wallis chi-squared = 17.171, df = 5, p-value = 0.004187

kw_pool_shannons <- kwfunction('Shannons', 'pool')
  ## Kruskal-Wallis chi-squared = 34.227, df = 5, p-value = 2.146e-06

kw_single_observed <- kwfunction('Observed', 'single')
  ## Kruskal-Wallis chi-squared = 28.164, df = 5, p-value = 3.381e-05

kw_single_shannons <- kwfunction('Shannons', 'single')
  ## Kruskal-Wallis chi-squared = 3.8534, df = 5, p-value = 0.5707

### Follow up with Dunn's test to see which pairwise tests are significant?
## Dunn test 


dunnfunction <- function(qfilt, batch){
  mangan_hill_df$Grouper <- paste(mangan_hill_df$CollectionMonth, mangan_hill_df$Site, sep="-")
  mangan_hill_df$Grouper <- as.factor(mangan_hill_df$Grouper)
  tmp <- dunnTest(Alpha_value ~ Grouper, 
                  data=mangan_hill_df %>% 
                    filter(BatchType == batch) %>% 
                    filter(Alpha_metric==qfilt), method = "bh") 
  tmp <- data.frame(tmp$res)
  tmp <- tmp %>% mutate(P.adj=round(P.adj, 3))
  tmp <- tmp %>% mutate(Z=round(Z, 3))
  tmp <- tmp %>% mutate(P.unadj=round(P.unadj, 3))
  tmp %>% arrange(P.adj)
}

dunn_single_Observed <- dunnfunction('Observed', 'single')
  ## almost all significant pairwise comparisons are between sites on June/July shoulders
  ## June generally had higher richness than September...
    ## sit pairs for BH-adjusted compas included:
      ## JuneEN:SeptHB, JuneHB:JuneEN, JuneEN:SeptHB, JuneHB:SeptHB, JuneEN:SeptEN, JulyEN:JulyHB

dunn_single_Shannons <- dunnfunction('Shannons', 'single')
  ## no pairwise significant values (as expected)

dunn_pool_Observed <- dunnfunction('Observed', 'pool')
  ## Just one pairwise comparison with p < 0.01... JulyHB:SeptEN. All other comps p ~>= 0.1

dunn_pool_Shannons <- dunnfunction('Shannons', 'pool')
  ## all significant pairwise comparisons are between two different sites
  ## unlike single samples, pooled samples have significant pairwise comparisons can vary between neighboring months (June:July, July:September)
  ## JulyHB:JulyEN, JuneHB:SeptEN, JulyHB:SeptEN, SeptEN:SeptHB, JulyHB:JuneEN all p <= 0.05 (BH-adjusted)



################################################################################
## unused code
################################################################################

# anovafunction <- function(qval, filename) {
#   aov.out <- aov(Alpha_value ~ CollectionMonth * Site, data=mangan_hill_df %>% filter(Alpha_metric == qval))
#   capture.output(summary(aov.out),file=paste0("~/Repos/mysosoup/data/text_tables/anovas/",filename))
# }
# anovafunction("Observed", "aov_mangan_alpha_q0.txt")
# anovafunction("Shannons", "aov_mangan_alpha_q1.txt")
# anovafunction("Simpsons", "aov_mangan_alpha_q2.txt")


# capture.output(kwfunction(mangan_hill_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q0.txt")
# capture.output(kwfunction(mangan_hill_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q1.txt")
# capture.output(kwfunction(mangan_hill_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q2.txt")

# capture.output(dunnfunction(mangan_hill_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q0.txt")
# capture.output(dunnfunction(mangan_hill_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q1.txt")
# capture.output(dunnfunction(mangan_hill_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q2.txt")
