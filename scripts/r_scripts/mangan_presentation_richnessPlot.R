library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)
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

## import taxa data, modify ambiguous_taxa as NA
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"
row.names(taxa) <- taxa$ASVid

## import metadata and add/modify labels for plot
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE) %>%
  select(-RoostDay, -RoostMonth, -RoostYear)
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$CollectionMonth[is.na(meta$CollectionMonth)] <- "control"
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "6", gsub("6", "June", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "7", gsub("7", "July", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "9", gsub("9", "September", meta$CollectionMonth), meta$CollectionMonth)
meta$ContamArea[is.na(meta$ContamArea)] <- "isolate"
meta$SiteMonth <- paste(meta$Site, meta$CollectionMonth, sep="-")
meta <- meta %>% rename(Month=CollectionMonth)

## import read data; filter out samples not included in Machine Learning analysis
## import dataset, filter out reads all non-Arthropod reads without (at least) Family-rank info, ..
# alternatively: download files:
## download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza", "tmp.qza")
# then set PATH to whatever directory needed:
# qzapath = "PATH/TO/tmp.qza"

## ..filter out remaining ASVs present exclusively in control samples, then calculate Hill Number values
## run from local: qzapath = "~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza"
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
data <- merge(df.tmp, meta)
rm(df.tmp, meta)

selectClass <- c("Arachnida", "Insecta")
order_by_month_Sumry <- data %>% 
  filter(class_name %in% selectClass) %>% 
  filter(!is.na(order_name)) %>%
  group_by(Site, Month, class_name, order_name) %>% 
  summarise(counts=n()) %>%
  group_by(Site, Month) %>%
  mutate(MonthlySum=sum(counts)) %>%
  mutate(FracMonthCounts=(counts/MonthlySum)*100)
order_less_than1 <- order_by_month_Sumry %>% filter(FracMonthCounts < 1)
remainder <- order_less_than1 %>% 
  group_by(Site, Month) %>% 
  summarise(counts=sum(counts)) %>%
  mutate(class_name="other", order_name="other")
remainder <- as.data.frame(remainder[,c("Site", "Month", "class_name", "order_name", "counts")])
order_more_than1 <- as.data.frame(order_by_month_Sumry %>% 
  filter(FracMonthCounts >= 1) %>%
  select(-FracMonthCounts, -MonthlySum))
plot_df <- rbind(remainder, order_more_than1)


pal9 <- c("#a14462", "#4c9568", "#f06152", "#76afda", "#ffc556", 
          "#eb998b", "#b0d45d", "#5066a1", "gray50")
plot_df$order_name <- factor(plot_df$order_name, levels = c("Araneae", "Coleoptera", "Diptera", 
                                                            "Hemiptera", "Hymenoptera", "Lepidoptera", 
                                                            "Trichoptera", "Psocodea", "other"))

plot_df$Month <- factor(plot_df$Month, 
                        levels = c("June", "July", "September"))

ggplot(plot_df %>% filter(!is.na(order_name)), aes(x=Month, y=counts, fill=order_name), color='black') +
  geom_bar(stat="identity") +
  facet_wrap(~ Site, nrow=2) + coord_flip() +
  scale_fill_manual(values=pal9) +
  labs(x="", y="observations of taxa", fill="Taxa Order") +
  theme_devon()


## summarize number of distinct ASVs per sample... distributions of ASVs per sample observed per Site/Month
ASVperSample <- data %>% group_by(Site, Month, SampleID) %>% summarise(nASVs=n_distinct(ASVid))

ASVperSample$Month <- factor(ASVperSample$Month, levels=c("June", "July", "September"))

ggplot(ASVperSample, aes(x=Month, y=nASVs)) + 
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(width=0.2, alpha=0.7) +
  facet_wrap(~ Site, nrow=2) +
  labs(x="", y="unique sequence variants per sample") +
  theme_devon()
