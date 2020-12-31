library(tidyverse)
library(vegan)
library(scales)
library(qiime2R)
library(reshape2)
library(formattable)
library(ggpubr)
library(Matrix)
library(multcompView)
library(FSA)
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
taxa <- read_csv(file="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/taxonomy/filtd_tax_dataframe_ALL.csv")
colnames(taxa)[1] <- "ASVid"

## import list of nonMYSO samples: 
nonMYSO_sampleID_list <- read_table(file="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/host/nonMYSO_sampleIDlist.txt", col_names = FALSE)

## import OTU table from QZA artifact
## download qza file:
# download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/Mangan.clust_p985_table_Rarefyd.qza", "Mangan.clust_p985_table_Rarefyd.qza")
# then set PATH to whatever directory needed:
# qzapath = "PATH/TO/tmp.qza"

## alternatively, run from local:
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
colnames(tmp) <- c("ASVid", "SampleID", "Reads")  ## note these are OTUs, not ASVs, but this matches the taxonomy info
df.tmp <- merge(tmp, taxa)
rm(tmp)
df.tmp <- merge(df.tmp, meta)
df.tmp <- df.tmp %>% filter(BatchType == "single" & !SampleID %in% nonMYSO_sampleID_list$X1)
mat.tmp <- dcast(data = df.tmp, formula = SampleID ~ ASVid, value.var='Reads', fill = 0)
row.names(mat.tmp) <- mat.tmp$SampleID
mat.tmp$SampleID <- NULL
tmp.hill <- data.frame(renyi(mat.tmp, scales = c(0,1), hill=TRUE)) %>% mutate(SampleID = row.names(.))
colnames(tmp.hill)[1:2] <- c("Observed", "Shannons")
Alpha_df <- gather(tmp.hill, key="Alpha_metric", value = "Alpha_value", c('Observed', 'Shannons'))
mangan_hill_df <- merge(Alpha_df, meta)
rm(tmp.hill, mat.tmp, Alpha_df, df.tmp, qzapath, nonMYSO_sampleID_list)


################################################################################
## Step 1 - test to see if the data is normal
## (hint, it's not)... going to use KW test for comparisons between Site-Month groups
## run test for each alpha metric (observed and Shannon's)
################################################################################

obs_alpha <- mangan_hill_df %>% 
  filter(Alpha_metric=="Observed") %>% 
  select(Alpha_value) %>% pull()

shan_alpha <- mangan_hill_df %>% 
  filter(Alpha_metric=="Shannons") %>% 
  select(Alpha_value) %>% pull()

shapiro.test(obs_alpha)
## significant result here indicates the data are VIOLATING assumption of normality. oops.
## W = 0.98457, p-value = 0.03592

shapiro.test(shan_alpha)
## W = 0.87179, p-value = 1.39e-11
## yeah, really not normal.


## see also the visualizations, which suggest an almost normal Observed dataset, but not for Shannons
ggdensity(obs_alpha, 
          main = "Density plot of Alpha val",
          xlab = "Observed OTU")
## pretty normal distribution?

ggdensity(shan_alpha, 
          main = "Density plot of Alpha val",
          xlab = "Observed OTU")
## very NOT normal distribution - negative binomial? chi squared? Disney-park ride?

ggqqplot(obs_alpha)
## nice. follows along expected

ggqqplot(shan_alpha)
## many samples with very low Shannons numbers, few outliers exceeding expected.


################################################################################
## Step 2a - compare ranked means between site-month groups with Kruskal-Wallis
## Step 2b - generate a figure that shows:
##   i. boxplot of Site-Month data, per Alpha_metric, as well as per-sample data (jittered points)
##  ii. Wilcoxon pairwise rank sum test showing differences between each Site-Month group
## iii. Kruskal-Wallis p-value at top of each plot
################################################################################

###### 2a. Run Kruskal Wallis:
kwfunction <- function(alpha_metric){
  mangan_hill_df$Grouper <- paste(mangan_hill_df$CollectionMonth, mangan_hill_df$Site, sep="-")
  mangan_hill_df$Grouper <- as.factor(mangan_hill_df$Grouper)
  kruskal.test(Alpha_value ~ Grouper, 
               data = mangan_hill_df %>% filter(Alpha_metric == alpha_metric))
}

kw_observed <- kwfunction('Observed')
kw_observed
## Kruskal-Wallis chi-squared = 25.389, df = 5, p-value = 0.0001172
## differences in ranks

kw_shannons <- kwfunction('Shannons')
kw_shannons
## Kruskal-Wallis chi-squared = 2.1737, df = 5, p-value = 0.8246
## no differences in ranks

###### 2b. Make boxplot but include significant pairwise comparisons
### this is a bit of a complicated plot, so we'll break it up into three parts

## 2bi. First part... split the dataset into Observed/Shannon's values
mangan_hill_df$CollectionMonth <- gsub("September", "Sept", mangan_hill_df$CollectionMonth)
mangan_hill_df$Grouper <- paste(mangan_hill_df$CollectionMonth, mangan_hill_df$Site, sep="-")

observed_dat <- mangan_hill_df %>% filter(Alpha_metric=="Observed")
observed_dat$Grouper <- factor(observed_dat$Grouper, levels = c(
  "June-EN", "June-HB", "July-EN", "July-HB", "Sept-EN", "Sept-HB"))

shannon_dat <- mangan_hill_df %>% filter(Alpha_metric=="Shannons")
shannon_dat$Grouper <- factor(shannon_dat$Grouper, levels = c(
  "June-EN", "June-HB", "July-EN", "July-HB", "Sept-EN", "Sept-HB"))

## 2bii. Second part... perform pairwise Wilcoxon rank sum tests (pairwise comparisons)
## use this function to generate lettered values differentiating every pairwise comp.
lpf <- function(data){
  attach(data)
  tmp_tri <- pairwise.wilcox.test(Alpha_value, Grouper, p.adjust.method = "BH")
  detach()
  
  h0_mat <- tmp_tri$p.value
  `June-EN`= c(1, NA, NA, NA, NA, NA)
  `Sept-HB` = c(NA, NA, NA, NA, 1)
  h0_mat <- cbind(h0_mat, `Sept-HB`)
  h0_mat <- rbind(`June-EN`, h0_mat)
  h0_mat <- Matrix::forceSymmetric(h0_mat,uplo="L")
  h0_mat <- as.matrix(h0_mat)
  h0_mat[is.na(h0_mat)] <- 1
  
  tmp_lmat <- multcompLetters(h0_mat, compare="<=", threshold=0.05, Letters=letters)
  tmp_labels <- data.frame(tmp_lmat$Letters) %>% 
    rename(Letters = tmp_lmat.Letters) %>% 
    mutate(Grouper = row.names(.))
}

observed_labels <- lpf(observed_dat)
shannon_labels <- lpf(shannon_dat)
  ## bother running for shannon because KW test was not significant?

plotdat_observed <- merge(observed_dat, observed_labels)
plotdat_shannons <- merge(shannon_dat, shannon_labels)


## use a modified function from above to retain the actual p-values for supplementary file:
wilcoxpvalfunction <- function(data, metricval){
  attach(data)
  tmp_tri <- pairwise.wilcox.test(Alpha_value, Grouper, p.adjust.method = "BH")
  detach()
  
  h0_mat <- data.frame(tmp_tri$p.value)
  h0_mat$Pair <- row.names(h0_mat)
  
  h0_mat %>% 
    pivot_longer(-Pair, names_to="Names", values_to = "BH-adjusted p-val") %>% 
    filter(!is.na(`BH-adjusted p-val`)) %>% 
    mutate(Names = sub("\\.", "-", Names),
           Metric = metricval,
           `BH-adjusted p-val` = round(`BH-adjusted p-val`, 3)) %>% 
    arrange(Metric, `BH-adjusted p-val`)
}

wilcox_data_obs <- wilcoxpvalfunction(observed_dat, "Observed")
wilcox_data_sha <- wilcoxpvalfunction(shannon_dat, "Shannons")

all_wilcoxon_pairwise <- rbind(wilcox_data_obs, wilcox_data_sha)
write_csv(all_wilcoxon_pairwise, 
          "~/github/mysosoup/data/text_tables/wilcoxon_pairwise/all_alpha_pairwise_Wilcoxon.csv")


## 2biii. Make individual plots for Observed/Shannon's, and add KW values from previous calculation in plot
## get the KW values for the annotation using this function:
getKWlabel_forplot_function <- function(kw_data){
  kw_degfree <- as.data.frame(kw_data$parameter) %>% pull() %>% as.character()
  kw_chistat <- as.data.frame(kw_data$statistic) %>% pull() %>% round(3) %>%  as.character()
  kw_pval <- as.data.frame(kw_data$p.value) %>% pull() %>% round(6) %>% as.character()
  paste0("Kruskal-Wallis   H(", kw_degfree,") = ", kw_chistat, ", p = ", kw_pval)
}

observed_KWlabel <- getKWlabel_forplot_function(kw_observed)
shannons_KWlabel <- getKWlabel_forplot_function(kw_shannons)

## now make each plot
## if using color, match scheme with beta diversity figure
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)

## observed plot...
## toggle on/off to make colored/not
observed_plot <- ggplot(plotdat_observed) +
  scale_color_manual(values = v3pal) +
  geom_boxplot(aes(x=Grouper, y=Alpha_value), outlier.shape = NA) +
  geom_jitter(aes(x=Grouper, y=Alpha_value, color=CollectionMonth, shape=Site), width=0.15, size=2) +    ## toggle on for color version
  #geom_jitter(aes(x=Grouper, y=Alpha_value, shape=Site), width=0.15, size=2) +   ## toggle on for black/white version
  facet_wrap(~ Alpha_metric) +
  labs(x="", y="Estimated diversity\n", color="Month") +
  geom_label(data = plotdat_observed,
             aes(x = Grouper, y = max(Alpha_value*1.1), label=Letters),
             colour = "white", fontface = "bold", fill="black", size=4) +
  annotate("text", x=2.5, y=75, label = observed_KWlabel, size=5) +
  theme_devon() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        #legend.position = "none") ## toggle on for black/white version
        legend.position="top")   ## toggle on for color version

shannon_plot <- ggplot(plotdat_shannons) +
  scale_color_manual(values = v3pal) +
  geom_boxplot(aes(x=Grouper, y=Alpha_value), outlier.shape = NA) +
  geom_jitter(aes(x=Grouper, y=Alpha_value, color=CollectionMonth, shape=Site), width=0.15, size=2) + ## toggle on for color version
  #geom_jitter(aes(x=Grouper, y=Alpha_value, shape=Site), width=0.15, size=2) +  ## toggle on for black/white version
  facet_wrap(~ Alpha_metric) +
  labs(x="", y="Estimated diversity\n", color="Month") +
  geom_label(data = plotdat_shannons,
             aes(x = Grouper, y = max(Alpha_value*1.1), label=Letters),
             colour = "white", fontface = "bold", fill="black", size=4) +
  annotate("text", x=2.5, y=35, label = shannons_KWlabel, size=5) +
  theme_devon() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.position="top")   ## toggle on for color version
#legend.position = "none") ## toggle on for black/white version

## stitch together the two plots and save
ggarrange(observed_plot, shannon_plot, nrow = 2, common.legend=TRUE)  ## toggle on for black/white version
ggsave("~github/mysosoup/figures/Figure_2_alphaDiv_col.png", width = 20, height = 20, units = "cm")
ggsave("~/github/mysosoup/figures/Figure_2_alphaDiv_col.svg", width = 20, height = 20, units = "cm")


#ggarrange(observed_plot, shannon_plot, nrow = 2)  ## toggle on for black/white version
#ggsave("~github/mysosoup/figures/Figure_2_alphaDiv_bw.png", width = 20, height = 20, units = "cm")
#ggsave("~/github/mysosoup/figures/Figure_2_alphaDiv_bw.svg", width = 20, height = 20, units = "cm")

################################################################################
################################################################################
######## unused code

## Dunn's test direct testing instead of Wilcoxon?
## note we're adjusting for multiple comparisons with Benjamini-Hochberg correction:

dunnfunction <- function(alpha_metric){
  mangan_hill_df$Grouper <- paste(mangan_hill_df$CollectionMonth, mangan_hill_df$Site, sep="-")
  mangan_hill_df$Grouper <- as.factor(mangan_hill_df$Grouper)
  tmp <- dunnTest(Alpha_value ~ Grouper, 
                  data=mangan_hill_df %>% 
                    filter(Alpha_metric==alpha_metric), method = "bh") 
  tmp <- data.frame(tmp$res)
  tmp <- tmp %>% mutate(P.adj=round(P.adj, 3))
  tmp <- tmp %>% mutate(Z=round(Z, 3))
  tmp <- tmp %>% mutate(P.unadj=round(P.unadj, 3))
  tmp %>% arrange(P.adj)
}

dunn_Observed <- dunnfunction('Observed')
## almost all significant pairwise comparisons are between sites on June/July shoulders
## June generally had higher richness than September...
## sit pairs for BH-adjusted compas included:
## JuneEN:SeptHB, JuneHB:JuneEN, JuneEN:SeptHB, JuneHB:SeptHB, JuneEN:SeptEN, JulyEN:JulyHB

dunn_Shannons <- dunnfunction('Shannons')
## no pairwise significant values (as expected)
