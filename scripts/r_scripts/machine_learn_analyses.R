library(qiime2R)
library(tidyverse)
library(ggrepel)
library(xkcdcolors)
library(reshape2)
library(scales)
library(gganimate)
library(ggpubr)

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

## see https://scikit-learn.org/stable/modules/ensemble.html#feature-importance for info on machine learning process
## import Feature Importance data, calculate cumulative sum, then join all data rames together
## first for rarefied data
rmonth.dat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/fi.Month.tsv", delim = "\t", col_names = TRUE)
rsite.dat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/fi.Site.tsv", delim = "\t", col_names = TRUE)
rsitemonth.dat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/fi.SiteMonth.tsv", delim = "\t", col_names = TRUE)
rbatch.dat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/fi.batch.tsv", delim = "\t", col_names = TRUE)

## next for non rarefied data
month.dat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/nonrare/nonrare_fi.month.tsv", delim = "\t", col_names = TRUE)
site.dat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/nonrare/nonrare_fi.site.tsv", delim = "\t", col_names = TRUE)
sitemonth.dat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/nonrare/nonrare_fi.sitemonth.tsv", delim = "\t", col_names = TRUE)
batch.dat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/nonrare/nonrare_fi.batch.tsv", delim = "\t", col_names = TRUE)


################################################################################
# 1. Create cumulative sum plots for rarefied and unrarefied data (and combined)
################################################################################


## calculate cummulative sums for each dataset and combine into single df
datprep.function <- function(data, LearnType, RareType) {
  colnames(data)[1] <- "ASVid"
  data$cumsum <- cumsum(data$importance)
  data %>% mutate(LearnType = LearnType) %>% mutate(RareType=RareType) %>% mutate(labelorder = row.names(.))
}

rbatch_df <- datprep.function(rbatch.dat, "batch", "rarefied")
rmonth_df <- datprep.function(rmonth.dat, "month", "rarefied")
rsite_df <- datprep.function(rsite.dat, "site", "rarefied")
rsitemonth_df <- datprep.function(rsitemonth.dat, "sitemonth", "rarefied")
batch_df <- datprep.function(batch.dat, "batch", "unrarefied")
month_df <- datprep.function(month.dat, "month", "unrarefied")
site_df <- datprep.function(site.dat, "site", "unrarefied")
sitemonth_df <- datprep.function(sitemonth.dat, "sitemonth", "unrarefied")

all_df <- rbind(rbatch_df, rmonth_df, rsite_df, rsitemonth_df, batch_df, month_df, site_df, sitemonth_df)
all_df$labelorder <- as.numeric(all_df$labelorder)
all_df$LearnType <- factor(all_df$LearnType, levels=c("month", "site", "sitemonth", "batch"))

##plot; export at 800x600; save as 'feature_importance_Rarefied.vs.Nonrarefied'
ggplot(all_df, aes(x=labelorder, y=cumsum)) + 
  geom_line() +
  theme_devon() +
  geom_hline(yintercept = 0.5, color="blue", size=0.5, linetype="dotted") +
  geom_hline(yintercept = 0.9, color="red", size=0.5, linetype="dotted") +
  facet_grid(LearnType ~ RareType) +
  labs(x="number of ASVs", y="fraction ") +
  theme_devon()


## plot just rarefied; save as 'feature_importance_RarefiedOnly'; export at 500x500
ggplot(data = all_df %>% filter(RareType == "rarefied"), aes(x=labelorder, y=cumsum)) + 
  geom_line() +
  theme_devon() +
  geom_hline(yintercept = 0.5, color="blue", size=0.5, linetype="dotted") +
  geom_hline(yintercept = 0.9, color="red", size=0.5, linetype="dotted") +
  facet_wrap(LearnType ~ ., ncol=2) +
  labs(x="number of ASVs", y="fraction ") +
  theme_devon()

## plot just unrarefied; save as 'feature_importance_UnarefiedOnly'; export at 500x500
ggplot(data = all_df %>% filter(RareType == "unrarefied"), aes(x=labelorder, y=cumsum)) + 
  geom_line() +
  theme_devon() +
  geom_hline(yintercept = 0.5, color="blue", size=0.5, linetype="dotted") +
  geom_hline(yintercept = 0.9, color="red", size=0.5, linetype="dotted") +
  facet_wrap(LearnType ~ ., ncol=2) +
  labs(x="number of ASVs", y="fraction ") +
  theme_devon()


################################################################################
# 2. Use read data to explore abundance and occurrence based trends at Order, Family, and ASV-levels based on most significant ASVs ..
# .. detected in the machine learning data above (per factor)
################################################################################

## import read datasets for rarefied and non rarefied data:
qzaimport.function <- function(qzapath){
  featuretable <- read_qza(qzapath)
  mat.tmp <- featuretable$data
  rm(featuretable)
  df.tmp <- as.data.frame(mat.tmp)
  rm(mat.tmp)
  df.tmp$OTUid <- rownames(df.tmp)
  rownames(df.tmp) <- NULL
  tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
  colnames(tmp) <- c("ASVid", "SampleID", "Reads")
  tmp
}

## download file: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/seqs/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza", "rarfyd.qza")
## download file: download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/seqs/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza", "nonfyd.qza")
# save path for each file.. for example: 
  #    rfypath="~/Downloads/rarfyd.qza"
  #    nonpath="~/Downloads/nonfyd.qza"

# alternative: run from local rfypath="/Users/do/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza"
# alternative: run from local nonpath="/Users/do/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza"

rarfyd_reads <- qzaimport.function(rfypath)
nonfyd_reads <- qzaimport.function(nonpath)

## add ASValias; because all rarefied ASVs are present in nonrarefied set, we'll use the ..
## order of the ASV by total read abundance to set the labels (so ASV-1 has the most total reads, ASV-2 the second most, etc.):
tmp1 <- nonfyd_reads %>% group_by(ASVid) %>%  summarise(nReads=sum(Reads)) %>% arrange(-nReads) %>% mutate(ASValias=paste0("ASV-", row.names(.))) %>% select(-nReads)
nonfyd_reads <- merge(nonfyd_reads, tmp1)
rarfyd_reads <- merge(rarfyd_reads, tmp1)

################################################################################
#### side note: the rarefied reads have slightly fewer ASVs than the non rarefied; these represent a tiny fraction overall
#### they also have slightly larger number of samples (285 to 277): see comparisons below...
## This comparison illustrates how negligible the additional ASVs present in the non rarefied samples are:
## get number of distinct ASVs and total number of reads per dataset
rarfyd_reads %>% summarise(uniqASVs=n_distinct(ASVid), totalReads=sum(Reads), uniqSamps=n_distinct(SampleID)) ## 2,573 unique ASVs; 1440400 total seqs; 277 samples
nonfyd_reads %>% summarise(uniqASVs=n_distinct(ASVid), totalReads=sum(Reads), uniqSamps=n_distinct(SampleID)) ## 2,657 unique ASVs; 7554794 total seqs; 285 samples
## which ASVs are private to non rarefied (which are dropped when rarefying)
uniq_nonrarfyd_asvs <- setdiff(nonfyd_reads$ASVid, rarfyd_reads$ASVid)  ## 84 different 
## how many reads do these 84 ASVs comprise in the entire dataset?
nonfyd_reads %>% filter(ASVid %in% uniq_nonrarfyd_asvs) %>% summarise(uniqs=n_distinct(ASVid), totalReads=sum(Reads))  ## 2272 total seqs (just 0.03%)
################################################################################


################################################################################
## the plots:
## 1. using the entire rarefied dataset, we plot the occurrence and abundances of all ASVs grouped at:
  #  Order 
  #  Family rank
## these are illustrated as heatmaps

## 2. using the rarefied dataset, we highlight select ASVs that provide the most predictive power at discrimnating by SiteWeek:
################################################################################

## For all projects, merge read data with taxa information and metadata:
## import taxa data
## add taxonomy information
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"
row.names(taxa) <- taxa$ASVid

## import metadata
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE) %>% 
  filter(SampleType == "sample") %>% select(SampleID, BatchType, Roost, CollectionMonth, Site)
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "6", gsub("6", "June", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "7", gsub("7", "July", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "9", gsub("9", "September", meta$CollectionMonth), meta$CollectionMonth)
meta$SiteMonth <- paste(meta$Site, meta$CollectionMonth, sep="-")

## merge the data
rfy_plotdat <- merge(rarfyd_reads, taxa) %>% merge(., meta)
## replace any NA in the species_name or genus_name with "Unassigned" - used in the slope plots later..
rfy_plotdat <- rfy_plotdat %>% 
  mutate_if(is.factor, as.character) %>% 
  mutate(genus_name = replace_na(genus_name, "unassigned")) %>% 
  mutate(Labeler = paste(ASValias, genus_name, sep=" | "))
  


################################################################################
## Gathering the Order-rank information: 
## Plotting both occurrence and abundance information

## 1. Order rank heatmap: all rarefied reads, all ASVs 
order_ocr_sumry <- rfy_plotdat %>% 
  group_by(order_name, SiteMonth) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>% 
  spread(SiteMonth, nSamples, fill = 0) %>% 
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nSamples") %>%
  group_by(SiteMonth) %>%
  mutate(mSamples=sum(nSamples)) %>% 
  mutate(pSamples=round((nSamples/mSamples),2))

order_abu_sumry <- rfy_plotdat %>% 
  group_by(order_name, SiteMonth) %>% 
  summarise(nReads = sum(Reads)) %>% 
  spread(SiteMonth, nReads, fill = 0) %>% 
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nReads") %>%
  group_by(SiteMonth) %>%
  mutate(mReads=sum(nReads)) %>% 
  mutate(pReads=round((nReads/mReads),2))

## group ocurrence and abundance information:
order_sumry <- merge(order_ocr_sumry, order_abu_sumry, by = c('order_name', 'SiteMonth')) %>% 
  mutate(splitter=SiteMonth) %>% separate(., col=splitter, into=c("Site", "Month"), sep="-")

## set levels for plot 
order_sumry$SiteMonth <- factor(order_sumry$SiteMonth, levels = c("EN-June", "EN-July", "EN-September", "HB-June", "HB-July", "HB-September"))
order_sumry$Month <- factor(order_sumry$Month, levels = c("June", "July", "September"))

## plot ALL occurrence data:
## save as 'heatmap_order_ocur_allOrders_rfydat', export at 800x800
hm_order_rfy_occur <- ggplot(order_sumry, aes(x=Month, y=order_name, fill=pSamples)) +
  geom_tile(color='gray70') + 
  #coord_equal() +
  facet_wrap(~ Site) +
  scale_x_discrete(labels = c("Jun", "Jul", "Sep")) +
  scale_fill_gradient(low = 'gray99', high = '#6f7c00', breaks=c(0,0.08,0.16)) +
  labs(x="", y="Arthropod Order\n", fill="fraction of\nSamples") +
  scale_y_discrete(position = "left") +
  theme_devon() +
  theme(legend.position = "top",
        axis.title = element_text(angle=0, hjust = .5, vjust=.5, face = "bold"))

## plot ALL abundance data:
## save as 'heatmap_order_abun_allOrders_rfydat', export at 800x800
hm_order_rfy_abund <- ggplot(order_sumry, aes(x=Month, y=order_name, fill=pReads)) +
  geom_tile(color='gray70') + 
  #coord_equal() +
  facet_wrap(~ Site) +
  scale_x_discrete(labels = c("Jun", "Jul", "Sep")) +
  scale_fill_gradient(low = 'gray99', high = '#7f5e00', breaks=c(0,0.3,0.6)) +
  labs(x="", y="\n", fill="fraction of\nSequences") +
  scale_y_discrete(position = "left") +
  theme_devon() +
  theme(legend.position = "top")

## can plot together, but keep the separate scales in view:
## save as "Order_sumry_heatmap_TwoScale"; export at 900x450
plot_grid(hm_order_rfy_occur, hm_order_rfy_abund, labels = c("A", "B"), ncol=2)

## plot together:
## reshape data for single plot to facet (doing this to keep constant "Fill" values in same legend)
allOrder_plotdat <- order_sumry %>% 
  select(Site, Month, order_name, pReads, pSamples) %>%
  melt(data = ., 
       id.vars = c("Month", "order_name", "Site"),
       value.name = "Fill_value",
       variable.name = "Value_type") %>% 
  mutate(Value_type = gsub("^p", "", .$Value_type)) %>%  
  mutate(FacetName = paste(Site, Value_type, sep = "-"))

## set levels
allOrder_plotdat$Month <- factor(order_sumry$Month, levels = c("June", "July", "September"))
allOrder_plotdat$FacetName <- factor(allOrder_plotdat$FacetName,
                                     levels = c('EN-Samples', 'HB-Samples',
                                                'EN-Reads', 'HB-Reads'))

## plot; save as 'Order_sumry_heatmap_oneScale'; export at 700x350
ggplot(allOrder_plotdat, 
       aes(x=Month, y=order_name, fill=Fill_value)) +
         geom_tile(color='gray70') +
         #coord_equal() +
         facet_grid(~ FacetName) +
         scale_x_discrete(labels = c("Jun", "Jul", "Sep")) +
         scale_fill_gradient(low = 'gray99', high = '#7f5e00') +
         labs(x="", y="Arthropod Order\n", fill="fraction of data") +
         scale_y_discrete(position = "left") +
         theme_devon()
       
## notice how different the impressions are between occurrence and abundance..
## ..diptera dominate abundance, but not occurrence

## going to merge Orders that have less than 1% of overall abundance across entire dataset into a single group, then replot
## summarize the per-Order number of reeds
keepOrders <-  order_sumry %>% 
  group_by(order_name) %>% 
  summarise(sum_nReads=sum(nReads), sum_mReads=sum(mReads)) %>% 
  mutate(p_totReads=sum_nReads/sum_mReads) %>% 
  filter(p_totReads >= 0.01) %>% 
  select(order_name) %>% mutate_if(is.factor, as.character) %>% pull()

## group all other taxa not in the `keepOrders' into single group`
tmp_order_others <- order_sumry %>% 
  filter(!order_name %in% keepOrders) %>% 
  group_by(SiteMonth) %>% 
  summarise(nReads=sum(nReads), mReads=sum(mReads)) %>% 
  mutate(pReads=nReads/mReads) %>%  
  mutate(order_name="others") %>% 
  mutate(nSamples=NA, mSamples=NA, pSamples=NA) %>% 
  mutate(splitter=SiteMonth) %>% separate(., col=splitter, into=c("Site", "Month"), sep="-")

## subset the desired Orders:
tmp_order_select <- order_sumry %>% filter(order_name %in% keepOrders)
## bring back together:
grouped_order_sumry <- rbind(tmp_order_others, tmp_order_select)

## set levels for y axis order:
grouped_order_sumry$order_name <- factor(grouped_order_sumry$order_name, levels = c(
  "Araneae", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Psocodea", "Trichoptera", "others"))

## same plot as above, but selecting just those more abundant Orders; save as 'heatmap_order_abun_selectOrders_rfydat', export at 600x600
ggplot(grouped_order_sumry, aes(x=Month, y=order_name, fill=pReads)) +
  geom_tile(color='gray70') + 
  coord_equal() +
  facet_wrap(~ Site) +
  scale_x_discrete(labels = c("Jun", "Jul", "Sep")) +
  scale_fill_gradient(low = 'gray99', high = '#7f5e00') +
  labs(x="", y="Taxa Order\n", fill="fraction of\nSequences \nper group") +
  scale_y_discrete(position = "left") +
  theme_devon()



################################################################################
## Performing the same assessment at the Family level gets tricky - over 200 families...
## Using Random Forest classifier (machine learning) to identify ASVs associated with distinct 'SiteMonth' groups 
## We'll focus on just the top (28) ASVs for this group, as these represent 50% of the cummulative sum of the "importance" in the model

## note we're going to then follow up with a series of plots that focus on the ASVs by Month (not SiteMonth)...
## .. as we'll notice that there just isn't that much variation for SiteMonth (similar trends, just subtle changes)
## One note: some of the changes depend on whether you look at abundance or occurrence data - so we'll show these in separate plots
################################################################################


## pulling select ASVs for "SiteMonth" machine learning output
selectASVs <- rsitemonth_df %>% filter(cumsum < 0.51) %>% select(ASVid) %>% pull()

## generate abundance summary data, filtering from SiteMonth list of potentially informative ASVs
ASV_abu_Sumry <- rfy_plotdat %>% filter(ASVid %in% selectASVs) %>%
  mutate(Taxa = paste(order_name, ASValias, sep = "-")) %>% 
  group_by(SiteMonth, Taxa, ASValias, Labeler) %>%
  summarise(nReads = sum(Reads)) %>%
  spread(SiteMonth, nReads, fill = 0) %>% 
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nReads") %>% 
  group_by(SiteMonth) %>%
  mutate(mReads=sum(nReads)) %>% 
  mutate(pReads=round((nReads/mReads),2)) %>% 
  mutate(Spliter=SiteMonth) %>% separate(data = ., col = Spliter, into=c("Site", "Month"), sep = "-")

## genearte occurrence information with same selected ASVs
ASV_ocr_sumry <- rfy_plotdat %>% filter(ASVid %in% selectASVs) %>%
  mutate(Taxa = paste(order_name, ASValias, sep = "-")) %>% 
  group_by(SiteMonth, Taxa, ASValias, Labeler) %>%
  summarise(nSamples=n_distinct(SampleID)) %>% 
  spread(SiteMonth, nSamples, fill = 0) %>% 
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nSamples") %>%
  group_by(SiteMonth) %>%
  mutate(mSamples=sum(nSamples)) %>% 
  mutate(pSamples=round((nSamples/mSamples),2))

## merge data sets togerther
ASV_sumry <- merge(ASV_abu_Sumry, ASV_ocr_sumry, by = c('Taxa', 'SiteMonth', 'ASValias', 'Labeler')) %>% 
  separate(., col = Taxa, into=c("order_name", "delete1", "delete2"), sep="-") %>% 
  select(-delete1, -delete2) %>% 
  mutate(Highlight="yes")

## just six orders remain; all were part of original 8 plotted in Order heatmap

## set levels for plot
ASV_sumry$Month <- factor(ASV_sumry$Month, levels = c("June", "July", "September"))
ASV_sumry$order_name <- factor(ASV_sumry$order_name, levels = c("Araneae","Coleoptera","Diptera","Hemiptera","Lepidoptera","Psocodea"))
#pal6 <- c('windows blue', 'macaroni and cheese', 'gray50', 'faded green', 'dusty purple', 'pale red')
pal6 <- c('#3778bf', '#efb435', 'gray50', '#7bb274', '#825f87', '#d9544d')

## plotting read abundances with GENUS + ASV as label; save as 'slopeplot_abund_rfydat_bySiteMonth', export at 800x1000
ml_sp_abu <- ggplot(data = ASV_sumry, aes(x = Month, y = pReads, group=ASValias, label = Labeler, color=order_name)) +
  facet_wrap( ~ Site, ncol=2) +
  geom_jitter(data = ASV_sumry %>% filter(order_name!="Diptera"), size=1, width=0.01) +
  geom_line(data = ASV_sumry %>% filter(order_name!="Diptera"), size=.9) +
  geom_jitter(data = ASV_sumry %>% filter(order_name=="Diptera"), size=1, width = 0.01) +
  geom_line(data = ASV_sumry %>% filter(order_name=="Diptera"), alpha=0.7) +
  scale_color_manual(values=pal6) +
  labs(x="", y="fraction of Sequences per ASV\n", color="Arthropod Order") +
  geom_label_repel(data = ASV_sumry %>% filter(Month=="June" & pReads > 0.02), 
                   aes(color=order_name), fill="white", nudge_x = -5, direction = "y", size=4, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  geom_label_repel(data = ASV_sumry %>% filter(Month=="September" & pReads > 0.02), 
                   aes(color=order_name), fill="white", nudge_x = 5, direction = "y", size=4, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  theme_devon() +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size=14),
        axis.text.y = element_text(size=14), axis.text.x = element_text(size=14), strip.text = element_text(size=16))

## plotting read ocurrences with GENUS + ASV as label; save as 'slopeplot_occur_rfydat_bySiteMonth', export at 1000x1000
ml_sp_ocr <- ggplot(data = ASV_sumry, aes(x = Month, y = pSamples, group=ASValias, label = Labeler, color=order_name)) +
  facet_wrap( ~ Site, ncol=2) +
  geom_jitter(data = ASV_sumry %>% filter(order_name!="Diptera"), size=1, width=0.01) +
  geom_line(data = ASV_sumry %>% filter(order_name!="Diptera"), size=.9) +
  geom_jitter(data = ASV_sumry %>% filter(order_name=="Diptera"), size=1, width = 0.01) +
  geom_line(data = ASV_sumry %>% filter(order_name=="Diptera"), alpha=0.7) +
  scale_color_manual(values=pal6) +
  labs(x="", y="fraction of Samples per ASV\n", color="Arthropod Order") +
  geom_label_repel(data = ASV_sumry %>% filter(Month=="June" & pSamples > 0.04), 
                   aes(color=order_name), fill="white", nudge_x = -5, direction = "y", size=4, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  geom_label_repel(data = ASV_sumry %>% filter(Month=="September" & pSamples > 0.04), 
                   aes(color=order_name), fill="white", nudge_x = 5, direction = "y", size=4, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  theme_devon() +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size=14),
        axis.text.y = element_text(size=14), axis.text.x = element_text(size=14), strip.text = element_text(size=16))

## plot all in a giant plot; save as "slopeplot_allDat_bySiteMonth"; export at 
ggarrange(ml_sp_ocr, ml_sp_abu, common.legend = TRUE, nrow = 2, labels = c("A", "B"))

## alternative to slopePlot: use animation
## just running a xy scatterplot with abundance/occurence info, with facet by Site and ...
## ..transition state is Month
sitmnth.ani <- ggplot(data = ASV_sumry, aes(y = pReads, x=pSamples, color=order_name, label = ASValias)) +
  geom_point() +
  facet_grid(~ Site) +
  scale_color_manual(values=pal6) +
  geom_label_repel(data = ASV_sumry %>% filter(pReads > 0 | pSamples > 1),
                   aes(color=order_name), fill="white", size = 2.25, force = 5, seed = 42,
                   segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  labs(y="fraction Reads", x="fraction Samples", color="") +
  theme_devon() +
  theme(legend.position = "top") +
  guides(color = guide_legend(nrow=1)) +
  transition_states(Month, transition_length = 1, state_length = 2) + ggtitle('Month: {closest_state}')

## render:
animate(sitmnth.ani, renderer = gifski_renderer(loop=TRUE), width=1000)
anim_save("~/Repos/mysosoup/figures/gifs/sitemonth_ASVs.gif")


################################################################################
## the changes really aren't that drastic for either occurrence or abundance data for most ASVs... this likely ..
## .. reflects the fact that most differnces are due to Month, not the combination of Site and Month
## Next plot does the same as above for ASV slopeplots, but pulls the ASVs from the Random Forest classifier specific for ..
## .. training on Month (not SiteMonth)
################################################################################


## pulling select ASVs for "Month" machine learning output
## note how we have selected the top 80% this time - there is more predictive power in these ASVs than in SiteMonth...
## we set the bar at 50% of cummulative importance for SiteMonth and got 28 ASVs.. there are 30 ASVs here, but they pull in 80% of cummulative importance
select_MonthASVs <- rmonth_df %>% filter(cumsum < 0.81) %>% select(ASVid) %>% pull()
## how many of these are the same ASVs? 
length(intersect(select_MonthASVs, selectASVs)) ## 22 are shared (among 28 and 30 for SiteMonth and Month, respectively)
  ## that's important - the two different models are converging on a similar set of winners

## generate abundance summary data, filtering from SiteMonth list of potentially informative ASVs
ASV_abu_Month_Sumry <- rfy_plotdat %>% filter(ASVid %in% select_MonthASVs) %>%
  mutate(Taxa = paste(order_name, ASValias, sep = "-")) %>% 
  group_by(CollectionMonth, Taxa, Labeler, ASValias, ASVid) %>%
  summarise(nReads = sum(Reads)) %>%
  spread(CollectionMonth, nReads, fill = 0) %>% 
  gather(June, July, September, key = "Month", value="nReads") %>% 
  group_by(Month) %>%
  mutate(mReads=sum(nReads)) %>% 
  mutate(pReads=round((nReads/mReads),2))

## genearte occurrence information with same selected ASVs
ASV_ocr_Month_Sumry <- rfy_plotdat %>% filter(ASVid %in% select_MonthASVs) %>%
  mutate(Taxa = paste(order_name, ASValias, sep = "-")) %>% 
  group_by(CollectionMonth, Taxa, Labeler, ASValias, ASVid) %>%
  summarise(nSamples=n_distinct(SampleID)) %>% 
  spread(CollectionMonth, nSamples, fill = 0) %>% 
  gather(June, July, September, key = "Month", value="nSamples") %>% 
  group_by(Month) %>%
  mutate(mSamples=sum(nSamples)) %>% 
  mutate(pSamples=round((nSamples/mSamples),2))

## merge data sets togerther
ASV_Month_sumry <- merge(ASV_abu_Month_Sumry, ASV_ocr_Month_Sumry, 
                         by=c('Taxa', 'Month', 'Labeler', 'ASVid', 'ASValias')) %>% 
  separate(., col = Taxa, into=c("order_name", "delete1", "delete2"), sep="-") %>% 
  select(-delete1, -delete2) %>% 
  mutate(Highlight="yes")

## set levels for plot
ASV_Month_sumry$Month <- factor(ASV_Month_sumry$Month, levels = c("June", "July", "September"))
ASV_Month_sumry$order_name <- factor(ASV_Month_sumry$order_name, levels = c("Araneae","Diptera","Hemiptera","Lepidoptera","Psocodea"))
  ## note we've lost the Coleopteran Order samples here... need to keep color scheme consistent

pal5 <- c('#3778bf', 'black', '#7bb274', '#825f87', '#d9544d')

## plotting read abundances with GENUS + ASV as label; save as 'slopeplot_abund_rfydat_byMonth', export at 1200x800
ggplot(data = ASV_Month_sumry, aes(x = Month, y = pReads, group=Labeler, label = Labeler, color=order_name)) +
  geom_point(data = ASV_Month_sumry %>% filter(order_name!="Diptera"), size=1, width=0.01) +
  geom_line(data = ASV_Month_sumry %>% filter(order_name!="Diptera"), size=.9) +
  geom_point(data = ASV_Month_sumry %>% filter(order_name=="Diptera"), size=1, width = 0.01) +
  geom_line(data = ASV_Month_sumry %>% filter(order_name=="Diptera"), alpha=0.7) +
  scale_color_manual(values=pal5) +
  labs(x="", y="fraction of Sequences per ASV\n", color="Taxa Order") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=14),
        axis.text.y = element_text(size=14), axis.text.x = element_text(size=14), strip.text = element_text(size=16)) +
  geom_label_repel(data = ASV_Month_sumry %>% filter(Month=="June" & pReads > 0.02), 
                   aes(color=order_name), fill="white", nudge_x = -5, direction = "y", size=4, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  geom_label_repel(data = ASV_Month_sumry %>% filter(Month=="September" & pReads > 0.02), 
                   aes(color=order_name), fill="white", nudge_x = 5, direction = "y", size=4, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  guides(color = guide_legend(nrow=1))


## plotting read ocurruences with GENUS + ASV as label; save as 'slopeplot_occur_rfydat_byMonth', export at 1200x800
ggplot(data = ASV_Month_sumry, aes(x = Month, y = pSamples, group=Labeler, label = Labeler, color=order_name)) +
  geom_point(data = ASV_Month_sumry %>% filter(order_name!="Diptera"), size=1, width=0.01) +
  geom_line(data = ASV_Month_sumry %>% filter(order_name!="Diptera"), size=.9) +
  geom_point(data = ASV_Month_sumry %>% filter(order_name=="Diptera"), size=.5, width = 0.01) +
  geom_line(data = ASV_Month_sumry %>% filter(order_name=="Diptera"), alpha=0.7) +
  scale_color_manual(values=pal5) +
  labs(x="", y="fraction of Samples per ASV\n", color="Taxa Order") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=14),
        axis.text.y = element_text(size=14), axis.text.x = element_text(size=14), strip.text = element_text(size=16)) +
  geom_label_repel(data = ASV_Month_sumry %>% filter(Month=="June" & pSamples > 0.02), 
                   aes(color=order_name), fill="white", nudge_x = -5, direction = "y", size=4, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  geom_label_repel(data = ASV_Month_sumry %>% filter(Month=="September" & pSamples > 0.02), 
                   aes(color=order_name), fill="white", nudge_x = 5, direction = "y", size=4, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  guides(color = guide_legend(nrow=1))


################################################################################
## one foray into gganimate to show the entire dataset moving by Month 
## plotting nReads ~ nSamples, with points as ASVs, and labels as ASValias 
require(gganimate)
# see: https://cran.r-project.org/web/packages/gganimate/vignettes/gganimate.html for more information

## generate per-ASV summary using all ASVs, then separate amongst MachineLearning taggeed and not...
## generate abundance summary data, filtering OUT SiteMonth list of potentially informative ASVs
ASV_abu_Month_Sumry_background <- rfy_plotdat %>% filter(!ASVid %in% select_MonthASVs) %>%
  mutate(Taxa = paste(order_name, ASValias, sep = "-")) %>% 
  group_by(CollectionMonth, Taxa, Labeler, ASValias, ASVid) %>%
  summarise(nReads = sum(Reads)) %>%
  spread(CollectionMonth, nReads, fill = 0) %>% 
  gather(June, July, September, key = "Month", value="nReads") %>% 
  group_by(Month) %>%
  mutate(mReads=sum(nReads)) %>% 
  mutate(pReads=round((nReads/mReads),2))

## genearte occurrence information with same selected ASVs
ASV_ocr_Month_Sumry_background <- rfy_plotdat %>% filter(!ASVid %in% select_MonthASVs) %>%
  mutate(Taxa = paste(order_name, ASValias, sep = "-")) %>% 
  group_by(CollectionMonth, Taxa, Labeler, ASValias, ASVid) %>%
  summarise(nSamples=n_distinct(SampleID)) %>% 
  spread(CollectionMonth, nSamples, fill = 0) %>% 
  gather(June, July, September, key = "Month", value="nSamples") %>% 
  group_by(Month) %>%
  mutate(mSamples=sum(nSamples)) %>% 
  mutate(pSamples=round((nSamples/mSamples),2))

## merge data together:
ASV_Month_Sumry_background <- merge(ASV_abu_Month_Sumry_background, ASV_ocr_Month_Sumry_background,
                                    by=c('Taxa', 'Month', 'Labeler', 'ASVid', 'ASValias')) %>% 
  separate(., col = Taxa, into=c("order_name", "delete1", "delete2"), sep="-") %>% 
  select(-delete1, -delete2) %>% 
  mutate(Highlight="no")

## combine with selected ASVs into single df
ASV_Month_sumry_fulldat <- rbind(ASV_Month_Sumry_background, ASV_Month_sumry)

## set levels for plot
ASV_Month_sumry_fulldat$Month <- factor(ASV_Month_sumry_fulldat$Month, levels = c("June", "July", "September"))

## plot the entire dataset, highlighting just the selected ASVs
full.ani <-ggplot(data = ASV_Month_sumry_fulldat, 
                  aes(x=nSamples, y=nReads, label=ASValias)) +
  geom_point(data = ASV_Month_sumry_fulldat %>% filter(Highlight=="no") %>% filter(order_name %in% nonselect_orderNames), 
             aes(color=order_name), size=3) +
  geom_point(data = ASV_Month_sumry_fulldat %>% filter(Highlight=="no"), color="gray90", size=1.5) +
  geom_point(data = ASV_Month_sumry_fulldat %>% filter(Highlight=="yes"), aes(color=order_name), size=4) +
  geom_label(data = ASV_Month_sumry_fulldat %>% filter(Highlight=="yes") %>% filter(nSamples > 20),
             aes(color=order_name), label.size = 0.25, size=5, nudge_x = -10) +
  scale_color_manual(values=pal5) +
  scale_y_continuous(labels=comma) +
  labs(x="Samples", y="Sequences", color="", 
       subtitle = "gray points reflect all ASVs in all Taxa not selected as significant to RF model") +
  guides(color = guide_legend(nrow=1)) +
  theme_devon() + 
  theme(legend.position = "top") +
  transition_states(Month, transition_length = 1, state_length = 4) + ggtitle('Month: {closest_state}')

## render:
animate(full.ani, renderer = gifski_renderer(loop=TRUE))
anim_save("~/Repos/mysosoup/figures/gifs/alldat_selectOrders_onepane.gif")


## now we plot just the top 8 Orders, grouping all other order_name as "others" 
## this let's us see how dynamic each order is separately, which we couldn't see in the other plot
# set palette; adding in hymenopteran and trichopteran to earlier palette schemes
pal8 <- c('#3778bf', '#efb435', 'black', '#7bb274', '#ff028d', '#825f87', '#d9544d', '#a87900')

## levels for plot
ASV_Month_sumry_fulldat$order_name <- factor(ASV_Month_sumry_fulldat$order_name, levels = c(
  "Araneae","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Psocodea","Trichoptera"))
ASV_Month_sumry_fulldat$Month <- factor(ASV_Month_sumry_fulldat$Month, levels = c("June", "July", "September"))

## ASValias label list
plotlist <- c('ASV-1', 'ASV-13', 'ASV-25', 'ASV-74', 'ASV-135', 'ASV-61', 'ASV-7', 'ASV-12', 'ASV-48', 'ASV-9', 'ASV-14', 'ASV-19')

## plot
byorder.ani <- ggplot(ASV_Month_sumry_fulldat %>% filter(order_name %in% keepOrders), 
                      aes(x=nSamples, y=nReads, color=order_name, label=Labeler)) +
  geom_point(data = ASV_Month_sumry_fulldat %>% filter(!ASVid %in% rmonth_df$ASVid) %>% filter(order_name %in% keepOrders),
             color="gray70", alpha=0.8) +
  geom_point(data = ASV_Month_sumry_fulldat %>% filter(ASVid %in% rmonth_df$ASVid) %>% filter(order_name %in% keepOrders), 
             aes(color=order_name), size=1.5, alpha=0.7) +
  geom_label(data = ASV_Month_sumry_fulldat %>% filter(ASValias %in% plotlist), nudge_x = -5, size = 4) +
  scale_color_manual(values=pal8) +
  facet_wrap(~ order_name, nrow=2) +
  scale_y_continuous(labels=comma, limits = c(0, 40000)) +
  labs(x="Samples", y="Sequences") +
  theme_devon() +
  theme(legend.position = "none") +
  guides(color = guide_legend(nrow=2)) +
  transition_states(Month, transition_length = 1, state_length = 4) + ggtitle('Month: {closest_state}')

animate(byorder.ani, renderer = gifski_renderer(loop=TRUE), width=1100)
anim_save("~/Repos/mysosoup/figures/gifs/alldat_selectOrders_byOrder.gif")


## final plot looks at just those crazy dipterans...
dip.ani <- ggplot(ASV_Month_sumry_fulldat %>% filter(order_name == "Diptera"), 
       aes(x=nSamples, y=nReads, label=Labeler)) +
  geom_point(data = ASV_Month_sumry_fulldat %>% filter(!ASVid %in% rmonth_df$ASVid) %>% filter(order_name == "Diptera"),
             color="gray70", alpha=0.8) +
  geom_point(data = ASV_Month_sumry_fulldat %>% filter(ASVid %in% rmonth_df$ASVid) %>% filter(order_name == "Diptera"), 
             color='black', size=1.5, alpha=0.7) +
  geom_label(data = ASV_Month_sumry_fulldat %>% filter(order_name == "Diptera") %>% filter(nSamples > 30 | nReads > 5000), 
             nudge_x = -2, size = 4) +
  scale_color_manual(values=pal8) +
  labs(x="Samples", y="Sequences") +
  theme_devon() +
  theme(legend.position = "none") +
  guides(color = guide_legend(nrow=2)) +
  transition_states(Month, transition_length = 1, state_length = 2) + ggtitle('Month: {closest_state}')
## render:
animate(dip.ani, renderer = gifski_renderer(loop=TRUE), width=1000)
anim_save("~/Repos/mysosoup/figures/gifs/dipteran_only.gif")


################################################################################
## unused plot info... unrarefied plot data produces plots for heatmaps just like rarefied
################################################################################

## merge the data
nonorfy_plotdat <- merge(nonfyd_reads, taxa) %>% merge(., meta)
## replace any NA in the species_name or genus_name with "Unassigned" - used in the slope plots later..
norfy_plotdat <- norfy_plotdat %>% 
  mutate_if(is.factor, as.character) %>% 
  mutate(genus_name = replace_na(genus_name, "unassigned")) %>% 
  mutate(Labeler = paste(ASValias, genus_name, sep=" | "))



################################################################################
## Gathering the Order-rank information: 
## Plotting both occurrence and abundance information

## 1. Order rank heatmap: all rarefied reads, all ASVs 
nofy_order_ocr_sumry <- norfy_plotdat %>% 
  group_by(order_name, SiteMonth) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>% 
  spread(SiteMonth, nSamples, fill = 0) %>% 
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nSamples") %>%
  group_by(SiteMonth) %>%
  mutate(mSamples=sum(nSamples)) %>% 
  mutate(pSamples=round((nSamples/mSamples),2))

nofy_order_abu_sumry <- norfy_plotdat %>% 
  group_by(order_name, SiteMonth) %>% 
  summarise(nReads = sum(Reads)) %>% 
  spread(SiteMonth, nReads, fill = 0) %>% 
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nReads") %>%
  group_by(SiteMonth) %>%
  mutate(mReads=sum(nReads)) %>% 
  mutate(pReads=round((nReads/mReads),2))

## group ocurrence and abundance information:
nofy_order_sumry <- merge(nofy_order_ocr_sumry, nofy_order_abu_sumry, by = c('order_name', 'SiteMonth')) %>% 
  mutate(splitter=SiteMonth) %>% separate(., col=splitter, into=c("Site", "Month"), sep="-")

## set levels for plot 
nofy_order_sumry$SiteMonth <- factor(nofy_order_sumry$SiteMonth, levels = c("EN-June", "EN-July", "EN-September", "HB-June", "HB-July", "HB-September"))
nofy_order_sumry$Month <- factor(nofy_order_sumry$Month, levels = c("June", "July", "September"))

## plot ALL occurrence data:
## save as 'heatmap_order_ocur_allOrders_rfydat', export at 800x800
nofy_hm_order_rfy_occur <- ggplot(nofy_order_sumry, aes(x=Month, y=order_name, fill=pSamples)) +
  geom_tile(color='gray70') + 
  #coord_equal() +
  facet_wrap(~ Site) +
  scale_x_discrete(labels = c("Jun", "Jul", "Sep")) +
  scale_fill_gradient(low = 'gray99', high = '#7f5e00') +
  labs(x="", y="Arthropod Order\n", fill="fraction of\nSamples") +
  scale_y_discrete(position = "right") +
  theme_devon() +
  theme(legend.position = "left")


## plot ALL abundance data:
## save as 'heatmap_order_abun_allOrders_rfydat', export at 800x800
nofy_hm_order_rfy_abund <- ggplot(nofy_order_sumry, aes(x=Month, y=order_name, fill=pReads)) +
  geom_tile(color='gray70') + 
  #coord_equal() +
  facet_wrap(~ Site) +
  scale_x_discrete(labels = c("Jun", "Jul", "Sep")) +
  scale_fill_gradient(low = 'gray99', high = '#7f5e00') +
  labs(x="", y=" \n", fill="fraction of\nSequences") +
  scale_y_discrete(position = "left") +
  theme_devon()

## can plot together, but keep the separate scales in view:
## save as "nonrarefied_Order_sumry_heatmap_TwoScale"; export at 1000x500
plot_grid(hm_order_rfy_occur, hm_order_rfy_abund, 
          labels = c("A", "B"))