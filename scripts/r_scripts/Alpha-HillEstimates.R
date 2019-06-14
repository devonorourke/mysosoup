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


## add taxonomy information
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"

## import dataset, filter out reads all non-Arthropod reads without (at least) Family-rank info, ..
# alternatively: download files:
## download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan.nonbatASVs.table.qza", "tmp.qza")
# then set PATH to whatever directory needed:
# qzapath = "PATH/TO/tmp.qza"

## ..filter out remaining ASVs present exclusively in control samples, then calculate Hill Number values
## run from local: qzapath = "~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza"
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
df.tmp <- merge(df.tmp, tinymeta)
onlyControlASVs <- setdiff(df.tmp %>% filter(SampleType=="control") %>% select(ASVid) %>% pull(), df.tmp %>% filter(SampleType=="sample") %>% select(ASVid) %>% pull())
df_filt.tmp <- df.tmp %>% filter(!ASVid %in% onlyControlASVs) %>% filter(phylum_name=="Arthropoda") %>% filter(!is.na(family_name)) %>% filter(SampleType == "sample")
rm(df.tmp)
mat.tmp <- dcast(data = df_filt.tmp, formula = SampleID ~ ASVid, value.var='Reads', fill = 0)
row.names(mat.tmp) <- mat.tmp$SampleID
mat.tmp$SampleID <- NULL
tmp.hill <- data.frame(renyi(mat.tmp, scales = c(0,1,2), hill=TRUE)) %>% mutate(SampleID = row.names(.))
colnames(tmp.hill)[1:3] <- c("q=0", "q=1", "q=2")
Alpha_df <- gather(tmp.hill, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
mangan_hill_df <- merge(Alpha_df, tinymeta)
rm(tmp.hill, mat.tmp, Alpha_df)


## plot parameter setup
mangan_hill_df$CollectionMonth <- factor(mangan_hill_df$CollectionMonth, levels = c("6", "7", "9"))
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)

## plot alpha diversity collectively
## save as all_Alpha_Hillvals; export at 900x697
ggplot(mangan_hill_df, aes(x=Labeler, y=Hill_value, color=CollectionMonth)) + 
  geom_jitter(width = 0.2, alpha=0.8) + 
  scale_color_manual(values = c(v3pal, "gray40"), labels=c("June", "July", "September")) +
  facet_grid( Hill_qType ~ .) +
  labs(x="", y="Estimated diversity", color = "Month") +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1),
        legend.position = "top")


##### Run ANOVA for group significance:
#### Run ANOVA for group significance:
anovafunction <- function(data, qval, filename) {
  aov.out <- aov(Hill_value ~ CollectionMonth * Site, data=data %>% filter(Hill_qType == qval))
  capture.output(summary(aov.out),file=paste0("~/Repos/mysosoup/data/text_tables/anovas/",filename))
}

anovafunction(mangan_hill_df, "q=0", "aov_mangan_alpha_q0.txt")
anovafunction(mangan_hill_df, "q=1", "aov_mangan_alpha_q1.txt")
anovafunction(mangan_hill_df, "q=2", "aov_mangan_alpha_q2.txt")


###### Run Kruskal Wallis for nonparametric test of same data:
kwfunction <- function(data, qval){
  data$Grouper <- paste(data$CollectionMonth, data$Site, sep="-")
  data$Grouper <- as.factor(data$Grouper)
  kruskal.test(Hill_value ~ Grouper, data = data %>% filter(Hill_qType == qval))
}

capture.output(kwfunction(mangan_hill_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q0.txt")
capture.output(kwfunction(mangan_hill_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q1.txt")
capture.output(kwfunction(mangan_hill_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q2.txt")


## Dunn test 
library(FSA)

dunnfunction <- function(data, qfilt){
  data$Grouper <- paste(data$CollectionMonth, data$Site, sep="-")
  data$Grouper <- as.factor(data$Grouper)
  dunnTest(Hill_value ~ Grouper, data=data %>% filter(Hill_qType==qfilt), method = "bh")
}

capture.output(dunnfunction(mangan_hill_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q0.txt")
capture.output(dunnfunction(mangan_hill_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q1.txt")
capture.output(dunnfunction(mangan_hill_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q2.txt")
