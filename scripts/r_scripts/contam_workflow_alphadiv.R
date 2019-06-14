library(tidyverse)
library(vegan)
library(scales)
library(qiime2R)
library(reshape2)

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
tinymeta <- meta %>% select(SampleID, Roost, CollectionMonth, Site, SampleType)
tinymeta$Site <- ifelse(tinymeta$Site == "Egner", gsub("Egner", "EN", tinymeta$Site), tinymeta$Site)
tinymeta$Site <- ifelse(tinymeta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", tinymeta$Site), tinymeta$Site)
tinymeta$CollectionMonth[is.na(tinymeta$CollectionMonth)] <- "control"
tinymeta$Labeler <- paste(tinymeta$Site, tinymeta$Roost, sep="-")
tinymeta$Labeler <- ifelse(tinymeta$Labeler == "control-control", gsub("control-control", "control", tinymeta$Labeler), tinymeta$Labeler)
tinymeta$CollectionMonth <- as.factor(tinymeta$CollectionMonth)

######################################################
# 1. Comparing Hill Numbers for all ASVs; no tax filtering
######################################################

## function will import dataset and calculate Hill Number values
## apply function to each dataset; set paths first.

alphaHill.function <- function(tablePath, filtType){
  featuretable <- read_qza(tablePath)
  mat.tmp <- featuretable$data
  t_mat.tmp <- t(mat.tmp)
  rm(featuretable, mat.tmp)
  tmp.hill <- data.frame(renyi(t_mat.tmp, scales = c(0,1,2), hill=TRUE)) %>% mutate(SampleID = row.names(.))
  colnames(tmp.hill)[1:3] <- c("q=0", "q=1", "q=2")
  Alpha_df <- gather(tmp.hill, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2')) %>% mutate(TableFiltType = filtType)
  merge(Alpha_df, tinymeta)
}

allASVpath = "~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.nonbatASVs.table.qza"
all_hill_df <- alphaHill.function(allASVpath, 'allASV')

## plot parameter setup
all_hill_df$CollectionMonth <- factor(all_hill_df$CollectionMonth, levels = c("6", "7", "9", "control"))
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)

## plotting all ASVs - notice no outliers in control samples - all highly abundant samples real samples 
## save as 'all_Alpha_Hillvals_allASVs'; export at 800x600
a <- ggplot(all_hill_df, aes(x=Labeler, y=Hill_value, color=CollectionMonth, label = SampleID)) + 
  geom_jitter(width = 0.2, alpha=0.8) + 
  scale_color_manual(values = c(v3pal, "gray40"), labels=c("June", "July", "September", "control")) +
  facet_wrap( ~ Hill_qType , scales = "free_y") +
  labs(x="", y="Estimated diversity", color = "Month") +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1),
        legend.position = "top")

## easier to see - coloring all samples read as those with 35 ASVs or less (same value is max number of ASVs in highest control sample)
ggplot(all_hill_df, aes(x=Labeler, y=Hill_value)) + 
  geom_jitter(data = all_hill_df %>% filter(Hill_value > 35), width = 0.2, alpha=0.5, color="black") + 
  geom_jitter(data = all_hill_df %>% filter(Hill_value <= 35), width = 0.2, alpha=0.5, color="red") + 
  facet_grid(Hill_qType ~ .) +
  labs(x="", y="Estimated diversity", color = "Month") +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1),
        legend.position = "top")

## Second way to view assess if ASVs are problematic: remove control samples, and discard any ASVs present in those samples from remaining dataset:
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
  merge(tmp, tinymeta)
}

dat <- tableimport.function(allASVpath)

## removing those ASVs in control samples from entire dataset
controlASVs <- dat %>% filter(SampleType == "control") %>% select(ASVid) %>% pull()
dat_filt <- dat %>% filter(!ASVid %in% controlASVs)
mat.tmp <- dcast(data = dat_filt, formula = SampleID ~ ASVid, value.var='Reads', fill = 0)
row.names(mat.tmp) <- mat.tmp$SampleID
mat.tmp$SampleID <- NULL
tmp.hill <- data.frame(renyi(mat.tmp, scales = c(0,1,2), hill=TRUE)) %>% mutate(SampleID = row.names(.))
colnames(tmp.hill)[1:3] <- c("q=0", "q=1", "q=2")
Alpha_df <- gather(tmp.hill, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
noNTCasv_df <- merge(Alpha_df, tinymeta)
rm(Alpha_df, tmp.hill, mat.tmp, dat_filt)

## plot this set of data
noNTCasv_df$CollectionMonth <- factor(noNTCasv_df$CollectionMonth, levels = c("6", "7", "9"))

## plotting remaining ASVs - notice the abundances match pretty much wdhat we saw before with all ASVs 
## save as 'all_Alpha_Hillvals_noNTCasvs'; export at 800x600
b <- ggplot(noNTCasv_df, aes(x=Labeler, y=Hill_value, color=CollectionMonth)) + 
  geom_jitter(width = 0.2, alpha=0.8) + 
  scale_color_manual(values = c(v3pal), labels=c("June", "July", "September")) +
  facet_grid(Hill_qType ~ .) +
  labs(x="", y="Estimated diversity", color = "Month") +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1),
        legend.position = "top")


## plot together:
## save as 'all_Alpha_Hillvals_allASVs_combo'; export at 1000x600
require(cowplot)
plot_grid(a, b, ncol=2)


######################################################
# 2. Comparing Hill Numbers for tax-filtered ASVs
######################################################

## creating dataset with control samples, but filtering out any ASV without ..
## .. at Arthropod classification containing at least Family-rank information 

## add taxonomy information
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"
df <- merge(dat, taxa)

## all samples
dat_taxfilt <- df %>% 
  filter(phylum_name == "Arthropoda") %>% 
  filter(!is.na(family_name))
## NTC-filtd ASVs 
taxfilt_NTCasvs <- dat_taxfilt %>% filter(SampleType == "control") %>% select(ASVid) %>% pull()
noNTCasvs_dat_filt <- dat_taxfilt %>% filter(!ASVid %in% taxfilt_NTCasvs)

## hill function 2.0 
hillfunction2 <- function(data){
  mat.tmp <- dcast(data = data, formula = SampleID ~ ASVid, value.var='Reads', fill = 0)
  row.names(mat.tmp) <- mat.tmp$SampleID
  mat.tmp$SampleID <- NULL
  tmp.hill <- data.frame(renyi(mat.tmp, scales = c(0,1,2), hill=TRUE)) %>% mutate(SampleID = row.names(.))
  colnames(tmp.hill)[1:3] <- c("q=0", "q=1", "q=2")
  Alpha_df <- gather(tmp.hill, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
  merge(Alpha_df, tinymeta)
}


hill_df_taxfilt_wNTCasv <- hillfunction2(dat_taxfilt)
hill_df_taxfilt_noNTCasv <- hillfunction2(noNTCasvs_dat_filt)

## plot tax-filtered data with NTC ASVs:
hill_df_taxfilt_wNTCasv$CollectionMonth <- factor(hill_df_taxfilt_wNTCasv$CollectionMonth, levels = c("6", "7", "9", "control"))

## save as 'taxfilt_Alpha_Hillvals_wNTCasvs'; export at 800x600
ggplot(hill_df_taxfilt_wNTCasv, aes(x=Labeler, y=Hill_value, color=CollectionMonth)) + 
  geom_jitter(width = 0.2, alpha=0.8) + 
  scale_color_manual(values = c(v3pal, 'gray40'), labels=c("June", "July", "September", "control")) +
  facet_wrap(~ Hill_qType, scales = "free_y") +
  labs(x="", y="Estimated diversity", color = "Month") +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1),
        legend.position = "top")


## plot tax filtered w/out NTC ASVs:
hill_df_taxfilt_noNTCasv$CollectionMonth <- factor(hill_df_taxfilt_noNTCasv$CollectionMonth, levels = c("6", "7", "9"))

## save as 'taxfilt_Alpha_Hillvals_noNTCasvs'; export at 800x600
ggplot(hill_df_taxfilt_noNTCasv, aes(x=Labeler, y=Hill_value, color=CollectionMonth)) + 
  geom_jitter(width = 0.2, alpha=0.8) + 
  scale_color_manual(values = c(v3pal, 'gray40'), labels=c("June", "July", "September", "control")) +
  facet_grid(Hill_qType ~ .) +
  labs(x="", y="Estimated diversity", color = "Month") +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1),
        legend.position = "top")


######################################################
# 3. Stat measures on all these alpha diversity values
######################################################

#### Run ANOVA for group significance:
anovafunction <- function(data, qval, filename) {
  aov.out <- aov(Hill_value ~ CollectionMonth * Site, data=data %>% filter(Hill_qType == qval))
  capture.output(summary(aov.out),file=paste0("~/Repos/mysosoup/data/text_tables/anovas/contam_checks/",filename))
}

## four input types; one anova per Hill Number (q = #):
anovafunction(all_hill_df, "q=0", "allASVs_wNTCasvs_q0.txt")
anovafunction(all_hill_df, "q=1", "allASVs_wNTCasvs_q1.txt")
anovafunction(all_hill_df, "q=2", "allASVs_wNTCasvs_q2.txt")

anovafunction(noNTCasv_df, "q=0", "allASVs_noNTCasvs_q0.txt")
anovafunction(noNTCasv_df, "q=1", "allASVs_noNTCasvs_q1.txt")
anovafunction(noNTCasv_df, "q=2", "allASVs_noNTCasvs_q2.txt")

anovafunction(hill_df_taxfilt_wNTCasv, "q=0", "taxfiltASVs_wNTCasvs_q0.txt")
anovafunction(hill_df_taxfilt_wNTCasv, "q=1", "taxfiltASVs_wNTCasvs_q1.txt")
anovafunction(hill_df_taxfilt_wNTCasv, "q=2", "taxfiltASVs_wNTCasvs_q2.txt")

anovafunction(hill_df_taxfilt_noNTCasv, "q=0", "taxfiltASVs_noNTCasvs_q0.txt")
anovafunction(hill_df_taxfilt_noNTCasv, "q=1", "taxfiltASVs_noNTCasvs_q1.txt")
anovafunction(hill_df_taxfilt_noNTCasv, "q=2", "taxfiltASVs_noNTCasvs_q2.txt")



###### Run Kruskal Wallis for nonparametric test of same data:
kwfunction <- function(data, qval){
  data$Grouper <- paste(data$CollectionMonth, data$Site, sep="-")
  data$Grouper <- as.factor(data$Grouper)
  kruskal.test(Hill_value ~ Grouper, data = data %>% filter(Hill_qType == qval))
}

capture.output(kwfunction(all_hill_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/allASVs_wNTCasvs_q0.txt")
capture.output(kwfunction(all_hill_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/allASVs_wNTCasvs_q1.txt")
capture.output(kwfunction(all_hill_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/allASVs_wNTCasvs_q2.txt")

capture.output(kwfunction(noNTCasv_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/allASVs_noNTCasvs_q0.txt")
capture.output(kwfunction(noNTCasv_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/allASVs_noNTCasvs_q1.txt")
capture.output(kwfunction(noNTCasv_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/allASVs_noNTCasvs_q2.txt")

capture.output(kwfunction(hill_df_taxfilt_wNTCasv, "q=0"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/taxfiltASVs_wNTCasvs_q0.txt")
capture.output(kwfunction(hill_df_taxfilt_wNTCasv, "q=1"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/taxfiltASVs_wNTCasvs_q1.txt")
capture.output(kwfunction(hill_df_taxfilt_wNTCasv, "q=2"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/taxfiltASVs_wNTCasvs_q2.txt")

capture.output(kwfunction(hill_df_taxfilt_noNTCasv, "q=0"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/taxfiltASVs_noNTCasvs_q0.txt")
capture.output(kwfunction(hill_df_taxfilt_noNTCasv, "q=1"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/taxfiltASVs_noNTCasvs_q1.txt")
capture.output(kwfunction(hill_df_taxfilt_noNTCasv, "q=2"),file="~/Repos/mysosoup/data/text_tables/kruskal/contam_checks/taxfiltASVs_noNTCasvs_q2.txt")


## Dunn test 
library(FSA)

dunnfunction <- function(data, qfilt){
  data$Grouper <- paste(data$CollectionMonth, data$Site, sep="-")
  data$Grouper <- as.factor(data$Grouper)
  dunnTest(Hill_value ~ Grouper, data=data %>% filter(Hill_qType==qfilt), method = "bh")
}

capture.output(dunnfunction(all_hill_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/allASVs_wNTCasvs_q0.txt")
capture.output(dunnfunction(all_hill_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/allASVs_wNTCasvs_q1.txt")
capture.output(dunnfunction(all_hill_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/allASVs_wNTCasvs_q2.txt")

capture.output(dunnfunction(noNTCasv_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/allASVs_noNTCasvs_q0.txt")
capture.output(dunnfunction(noNTCasv_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/allASVs_noNTCasvs_q1.txt")
capture.output(dunnfunction(noNTCasv_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/allASVs_noNTCasvs_q2.txt")

capture.output(dunnfunction(hill_df_taxfilt_wNTCasv, "q=0"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/taxfiltASVs_wNTCasvs_q0.txt")
capture.output(dunnfunction(hill_df_taxfilt_wNTCasv, "q=1"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/taxfiltASVs_wNTCasvs_q1.txt")
capture.output(dunnfunction(hill_df_taxfilt_wNTCasv, "q=2"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/taxfiltASVs_wNTCasvs_q2.txt")

capture.output(dunnfunction(hill_df_taxfilt_noNTCasv, "q=0"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/taxfiltASVs_noNTCasvs_q0.txt")
capture.output(dunnfunction(hill_df_taxfilt_noNTCasv, "q=1"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/taxfiltASVs_noNTCasvs_q1.txt")
capture.output(dunnfunction(hill_df_taxfilt_noNTCasv, "q=2"),file="~/Repos/mysosoup/data/text_tables/dunn/contam_check/taxfiltASVs_noNTCasvs_q2.txt")
