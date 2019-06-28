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
## We'll focus on just the top ASVs for this group, represented by those features in the 75th percentile of Importance

## note we're going to then follow up with a series of plots that focus on the ASVs by Month (not SiteMonth)...
## .. as we'll notice that there just isn't that much variation for SiteMonth (similar trends, just subtle changes)
## One note: some of the changes depend on whether you look at abundance or occurrence data - so we'll show these in separate plots
################################################################################


## pulling select ASVs for "SiteMonth" machine learning output
## selecting features in 75th percentile of "Importance" 
filt75_sm <- quantile(rsitemonth_df$importance, .75)
selectASVs <-rsitemonth_df %>% filter(importance >= filt75_sm) %>% pull(ASVid)

## generate abundance summary data for all data (no ASV filtering yet)
ASV_abu_Sumry <- rfy_plotdat %>% 
  mutate(Taxa = paste(order_name, ASValias, sep = "-")) %>% 
  group_by(SiteMonth, Taxa, ASValias, ASVid, Labeler) %>%
  summarise(nReads = sum(Reads)) %>%
  spread(SiteMonth, nReads, fill = 0) %>%
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nReads") %>%
  group_by(SiteMonth) %>%
  mutate(mReads=sum(nReads)) %>%
  mutate(pReads=round((nReads/mReads),2)) %>%
  mutate(Spliter=SiteMonth) %>% separate(data = ., col = Spliter, into=c("Site", "Month"), sep = "-")

## genearte occurrence information with same selected ASVs (no ASV filtering yet either)
nSampDat <- rfy_plotdat %>%
  select(SampleID, SiteMonth) %>% 
  group_by(SiteMonth) %>% 
  summarise(smSamples=n_distinct(SampleID))

ASV_ocr_sumry <- rfy_plotdat %>%
  mutate(Taxa = paste(order_name, ASValias, sep = "-")) %>%
  group_by(SiteMonth, Taxa, ASValias, Labeler) %>%
  summarise(nSamples=n_distinct(SampleID)) %>%
  spread(SiteMonth, nSamples, fill = 0) %>%
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nSamples") %>% 
  merge(., nSampDat) %>%
  mutate(pSamples=round((nSamples/smSamples),2))


## merge data sets togerther (still no ASV filtering)
ASV_sumry <- merge(ASV_abu_Sumry, ASV_ocr_sumry, by = c('Taxa', 'SiteMonth', 'ASValias', 'Labeler')) %>% 
  separate(., col = Taxa, into=c("order_name", "delete1", "delete2"), sep="-") %>% 
  select(-delete1, -delete2) %>% 
  mutate(Highlight="yes")

## now create filtered dataset
ASV_sumry_SLfiltd <- ASV_sumry %>% filter(ASVid %in% selectASVs)

## which orders would remain if we applied supervised learning filter?
OrderKeepers <- ASV_sumry_SLfiltd %>% distinct(order_name) %>% pull(order_name)
## eight orders would remain

## set levels for plot
ASV_sumry_SLfiltd$Month <- factor(ASV_sumry_SLfiltd$Month, levels = c("June", "July", "September"))
ASV_sumry_SLfiltd$order_name <- factor(ASV_sumry_SLfiltd$order_name, 
                                       levels = c("Araneae","Coleoptera","Diptera","Ephemeroptera", "Hemiptera","Lepidoptera","Psocodea","Trichoptera"))
pal8 <- c('#3778bf', '#efb435', 'black', '#7bb274', '#ff028d', '#825f87', '#d9544d', '#a87900')

## plotting read abundances with GENUS + ASV as label; save as 'slopeplot_abund_rfydat_bySiteMonth', export at 1000x1000
ml_sp_abu <- ggplot(data = ASV_sumry_SLfiltd,
                    aes(x = Month, y = pReads, group=ASValias, label = Labeler, color=order_name)) +
  facet_wrap( ~ Site, ncol=2) +
  geom_point(data = ASV_sumry_SLfiltd %>% filter(order_name!="Diptera"), size=1) +
  geom_line(data = ASV_sumry_SLfiltd %>% filter(order_name!="Diptera"), size=.9) +
  geom_point(data = ASV_sumry_SLfiltd %>% filter(order_name=="Diptera"), size=1) +
  geom_line(data = ASV_sumry_SLfiltd %>% filter(order_name=="Diptera"), alpha=0.7) +
  scale_color_manual(values=pal8) +
  scale_x_discrete(expand = c(0,2), labels = c("June", "July", "Sept.")) +
  labs(x="", y="fraction of Sequences per ASV per group\n", color="Arthropod Order") +
  geom_label_repel(data = ASV_sumry_SLfiltd %>% filter(Month=="June" & pReads > 0.02), 
                   aes(color=order_name), fill="white", nudge_x = -5, direction = "y", size=3, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  geom_label_repel(data = ASV_sumry_SLfiltd %>% filter(Month=="September" & pReads > 0.02), 
                   aes(color=order_name), fill="white", nudge_x = 5, direction = "y", size=3, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  theme_devon() +
  guides(color = guide_legend(nrow = 2)) +
  theme(legend.position = "top")

## plotting read ocurrences with GENUS + ASV as label; save as 'slopeplot_occur_rfydat_bySiteMonth', export at 1000x1000
ml_sp_ocr <- ggplot(data = ASV_sumry_SLfiltd, 
                    aes(x = Month, y = pSamples, group=ASValias, label = Labeler, color=order_name)) +
  geom_point() +
  facet_wrap( ~ Site, ncol=2) +
  geom_point(data = ASV_sumry_SLfiltd %>% filter(order_name!="Diptera"), size=1) +
  geom_line(data = ASV_sumry_SLfiltd %>% filter(order_name!="Diptera"), size=.9) +
  geom_point(data = ASV_sumry_SLfiltd %>% filter(order_name=="Diptera"), size=1) +
  geom_line(data = ASV_sumry_SLfiltd %>% filter(order_name=="Diptera"), alpha=0.7) +
  scale_color_manual(values=pal8) +
  scale_x_discrete(expand = c(0,2), labels = c("June", "July", "Sept.")) +
  labs(x="", y="fraction of Samples per ASV per group\n", color="Arthropod Order") +
  geom_label_repel(data = ASV_sumry_SLfiltd %>% filter(Month=="June" & pSamples > 0.20), 
                   aes(color=order_name), fill="white", nudge_x = -5, direction = "y", size=3, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  geom_label_repel(data = ASV_sumry_SLfiltd %>% filter(Month=="September" & pSamples > 0.20), 
                   aes(color=order_name), fill="white", nudge_x = 5, direction = "y", size=3, segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  theme_devon() +
  guides(colour = guide_legend(nrow = 2)) +
  theme(legend.position = "top")

## plot all in a giant plot; save as "slopeplot_allDat_bySiteMonth"; export at 1800x900
ggarrange(ml_sp_ocr, ml_sp_abu, common.legend = TRUE, labels = c("A", "B"))

## alternative to slopePlot to grasp more data at once: to Site/Month into unique facets and ..
## set sample detections and fractions of read abundances as X/Y scatterplot
## save as 'ml_pAbu_by_pOcr_scatterplot_byMonth_andSite'; export at 1000x1000

## split up the dataset into three parts:
## 1. ML data (solid circles, colored by 8 orders)
## 2. non ML data in same 8 orders as #1
## 3. non ML data in orders other than #2

MLasvs <- ASV_sumry_SLfiltd$ASValias %>% unique(.)
nonML_selectOrders_ASVs <- ASV_sumry %>% filter(order_name %in% OrderKeepers) %>% filter(!ASValias %in% ASV_sumry_SLfiltd$ASValias) %>% pull(ASValias) %>% unique(.)
nonML_otherOrders_ASVs <- ASV_sumry %>% filter(!order_name %in% OrderKeepers) %>% pull(ASValias) %>% unique(.)

## use ASV_sumry data.frame but rename "order_name" to "plotOrder_name" and subsitute all other Orders not in OrderKeeper to "Other"
ASV_sumry$plotOrder_name <- ASV_sumry$order_name
x <- unique(ASV_sumry$order_name) ## all Order names
y <- setdiff(x, OrderKeepers) ## all "Other" Order names
pats <- c('Blattodea|Entomobryomorpha|Hymenoptera|Isopoda|Mecoptera|Mesostigmata|Neuroptera|Odonata|Orthoptera|Poduromorpha|Sarcoptiformes|Strepsiptera|Trombidiformes')
ASV_sumry$plotOrder_name <- str_replace_all(ASV_sumry$plotOrder_name, pats, "Other")  ## this does the substitution

## use ASV_sumry data.frame, but add in field that distinguishes ML-labeled ASVs or not w/ TRUE/FALSE statement
ASV_sumry <- ASV_sumry %>% mutate(MLsample = ASValias %in% MLasvs)
## add in filterer for later faceting probs
ASV_sumry$MLsample_filter <- ASV_sumry$MLsample
ASV_sumry$MLsample_filter <- ifelse(ASV_sumry$MLsample_filter==TRUE, "SL-feature", "not-SL-feature")

## add in new color to 8-color palette to reflect "Other" now...
pal9 <- c(pal8, "gray60")

## set levels so that "Other" is placed at bottom of legend:
ASV_sumry$plotOrder_name <- factor(ASV_sumry$plotOrder_name, levels = c(
  "Araneae", "Coleoptera","Diptera","Ephemeroptera","Hemiptera","Lepidoptera","Psocodea","Trichoptera","Other"))
ASV_sumry$Month <- factor(ASV_sumry$Month, levels=c("June", "July", "September"))


## plot; save as 'ml_pAbu_by_pOcr_scatterplot_byMonth_andSite'; export at 900x900
## removed labels to avoid artibrary text labeling (we're already discriminating among ML and non-ML sites...)
ggplot() +
  geom_point(data = ASV_sumry,
             aes(y = pReads, x=pSamples, color=plotOrder_name, shape=MLsample),
             alpha = 0.8, size=2.5) +
  scale_shape_manual(values=c(0,19)) +
  facet_grid(Month ~ Site) +
  scale_color_manual(values=pal9) +
  scale_x_continuous(expand = c(0,0.15)) +
  geom_text_repel(data = ASV_sumry %>% filter(pSamples > 0.37),
                  aes(y = pReads, x=pSamples, color=order_name, label=Labeler), 
                  nudge_x = 5, force=5, direction="y", segment.size = 0.2, segment.alpha = 0.5, seed = 42, size=3) +
  #geom_label_repel(data = ASV_sumry %>% filter(pSamples > 0.4),
  #                 aes(y = pReads, x=pSamples, color=order_name, label=Labeler), 
  #                 fill="white", size = 2, force = 5, seed = 42, segment.size = 0.2, segment.alpha=0.5, direction="y", nudge_x = 10) +
  labs(y="fraction Reads", x="fraction Samples", color="Arthropod Order", shape="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=12)) +
  guides(color = guide_legend(nrow=2),
         shape=FALSE)

## can make these facets easier to identify differences with animation
## same plot as above, but eliminating Month facet with animation
## ..transition state is Month
sitmnth.ani <- ggplot() +
  geom_point(data = ASV_sumry,
             aes(y = pReads, x=pSamples, color=plotOrder_name, shape=MLsample),
             alpha = 0.8, size=2.5) +
  scale_shape_manual(values=c(0,19)) +
  facet_grid( ~ Site) +
  scale_color_manual(values=pal9) +
  scale_x_continuous(expand = c(0.05,0.15)) +
  geom_text_repel(data = ASV_sumry %>% filter(MLsample==TRUE) %>% filter(pReads > 0 | pSamples > 1),
                  aes(y = pReads, x=pSamples, color=order_name, label=Labeler), size = 3, force = 7, seed = 42,
                  segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  labs(y="fraction Reads", x="fraction Samples", color="Arthropod Order", shape="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=12), plot.title = element_text(size=22)) +
  guides(color = guide_legend(nrow=2), shape=FALSE) +
  transition_states(Month, transition_length = 1, state_length = 2) + ggtitle('Month: {closest_state}')

## render:
animate(sitmnth.ani, renderer = gifski_renderer(loop=TRUE), width=1000)
anim_save("~/Repos/mysosoup/figures/gifs/sitemonth_ASVs.gif")

##########
## animation shows lots of changes, but largely confined to Dipterans. Want to make a faceted plot looking at just Dipteran changes, as a function of ..
## .. sequence abundance and observations
## plot; save as ml_pAbu_by_pOcr_scatterplot_byMonth_andSite_DipteraOnly; export at 900x900

## set levels
ASV_sumry$MLsample_filter <- factor(ASV_sumry$MLsample_filter, levels = c("SL-feature", "not-SL-feature"))
ASV_sumry$Month <- factor(ASV_sumry$Month, levels = c("June", "July", "September"))

ggplot() +
  geom_point(data = ASV_sumry %>% filter(order_name=="Diptera"),
             aes(y = pReads, x=pSamples, shape=MLsample),
             alpha = 0.8, size=2.5, color="black") +
  scale_shape_manual(values=c(0,19)) +
  facet_grid(Month ~ Site) +
  scale_x_continuous(expand = expand_scale(add = 0.36)) +
  geom_label_repel(data = ASV_sumry %>% filter(order_name == "Diptera") %>% filter(MLsample==TRUE) %>% filter(pSamples > 0.3),
                   aes(y = pReads, x=pSamples, label=Labeler),
                   #color='black', fill="white", size = 2, force = 5, seed = 42,
                   color='black', fill="white", size = 2.25, force = 5, seed = 42, nudge_x = 5, direction = "y",
                   segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  geom_label_repel(data = ASV_sumry %>% filter(order_name == "Diptera") %>% filter(MLsample==FALSE) %>% filter(pReads > 0),
                   aes(y = pReads, x=pSamples, label=Labeler),
                   color='black', fill="white", size = 2.25, force = 5, seed = 42, nudge_x = -5, direction = "y",
                   segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  labs(y="fraction Reads", x="fraction Samples", color="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=16),
        plot.title = element_text(size=20, face = "bold")) +
  guides(color = FALSE, shape=FALSE)

## Animation for Dipteran only:
dipAni <- ggplot() +
  geom_point(data = ASV_sumry %>% filter(order_name=="Diptera"),
             aes(y = pReads, x=pSamples, shape=MLsample),
             alpha = 0.8, size=3, color="black") +
  scale_shape_manual(values=c(0,19)) +
  facet_grid( ~ Site) +
  scale_x_continuous(expand = expand_scale(add = 0.36)) +
  geom_label_repel(data = ASV_sumry %>% filter(order_name == "Diptera") %>% filter(MLsample_filter=='SL-feature') %>% filter(pSamples > 0.3),
                   aes(y = pReads, x=pSamples, label=Labeler),
                   color='black', fill="white", size = 4.5, force = 5, seed = 42, nudge_x = 5, direction="y",
                   segment.size = 0.2, segment.alpha=0.5, segment.colour = "gray50") +
  geom_label_repel(data = ASV_sumry %>% filter(order_name == "Diptera") %>% filter(MLsample_filter=='not-SL-feature') %>% filter(pReads > 0),
                   aes(y = pReads, x=pSamples, label=Labeler),
                   color='black', fill=NA, size = 4.5, force = 5, seed = 42,nudge_x = -5, direction="y",
                   segment.size = 0.3, segment.alpha=0.5, segment.colour = "black") +
  labs(y="fraction Reads", x="fraction Samples", color="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=16), plot.title = element_text(size=22, face = "bold")) +
  guides(color = FALSE, shape=FALSE) +
  transition_states(Month, transition_length = 1, state_length = 2) + ggtitle('Month: {closest_state}')

## render:
animate(dipAni, renderer = gifski_renderer(loop=TRUE), width=1000, height=500)
anim_save("~/Repos/mysosoup/figures/gifs/sitemonth_ASVs_DipteraOnly.gif")


##########
## alternate animation will focus on the changes within each Order at each Site for each Month
## removed most labels - too little
## selected single ASVs per Order
ASVs2labelbc <- c('ASV-1', 'ASV-25', 'ASV-3', 'ASV-60','ASV-7', 'ASV-9', 'ASV-14','ASV-19','ASV-81')

sitmnt_byOrder.ani <- ggplot() +
  scale_shape_manual(values=c(0,19)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  geom_label_repel(data = ASV_sumry %>% filter(ASValias %in% ASVs2labelbc),
                   aes(y = pReads, x=pSamples, color=plotOrder_name, label=Labeler), size = 2.75, seed = 42, nudge_y=0.03,
                   segment.size = 0.4, segment.alpha=0.9, segment.colour = "black", fill=NA) +
  geom_point(data = ASV_sumry,
             aes(y = pReads, x=pSamples, color=plotOrder_name, shape=MLsample),
             alpha = 0.8, size=2.5) +
  facet_grid(Site ~ plotOrder_name) +
  scale_color_manual(values=pal9) +
  labs(y="fraction Reads", x="fraction Samples", color="Arthropod Order", shape="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=12), plot.title = element_text(size=22)) +
  guides(color = guide_legend(nrow=2), shape=FALSE) +
  transition_states(Month, transition_length = 1, state_length = 2) + ggtitle('Month: {closest_state}')

## render:
animate(sitmnt_byOrder.ani, renderer = gifski_renderer(loop=TRUE), width=1000, height=500)
anim_save("~/Repos/mysosoup/figures/gifs/sitemonth_byOrder_ASVs.gif")



#########
## how many unique species are there in those Families targeted by SL model?
rfy_plotdat %>% 
  filter(ASVid %in% selectASVs) %>%  ## drop this if you don't want to focus on 75th percentile
  group_by(order_name, family_name) %>%
  summarise(nTaxa=n_distinct(nTaxa=ASValias)) %>%
  mutate(pTaxa=round(nTaxa/(sum(nTaxa)),2)) %>% 
  arrange(-nTaxa)

## how many species in Lep family, and top 4 Dips?
rfy_plotdat %>% filter(family_name == "Tortricidae" & order_name == "Lepidoptera") %>% distinct(species_name)
## 36 + 1 species named
rfy_plotdat %>% filter(family_name == "Culicidae" & order_name == "Diptera") %>% distinct(species_name)
## 16 species named
rfy_plotdat %>% filter(family_name == "Limoniidae" & order_name == "Diptera") %>% distinct(species_name)
## 12 species named
rfy_plotdat %>% filter(family_name == "Tipulidae" & order_name == "Diptera") %>% distinct(species_name)
## 12 species named
rfy_plotdat %>% filter(family_name == "Chironomidae" & order_name == "Diptera") %>% distinct(genus_name)
## 18 species named 

## how many ASVs?
rfy_plotdat %>% filter(family_name == "Tortricidae" & order_name == "Lepidoptera") %>% distinct(ASValias)
## 119 ASVs
rfy_plotdat %>% filter(family_name == "Culicidae" & order_name == "Diptera") %>% distinct(ASValias)
## 255 ASVs
rfy_plotdat %>% filter(family_name == "Limoniidae" & order_name == "Diptera") %>% distinct(ASValias)
## 252 ASVs
rfy_plotdat %>% filter(family_name == "Tipulidae" & order_name == "Diptera") %>% distinct(ASValias)
## 137 ASVs
rfy_plotdat %>% filter(family_name == "Chironomidae" & order_name == "Diptera") %>% distinct(ASValias)
## 244 ASVs


## which ASVs are private to each Site?
ENasvs <- rfy_plotdat %>% filter(Site=="EN") %>% select(Labeler) %>% pull()
HBasvs <- rfy_plotdat %>% filter(Site=="HB") %>% select(Labeler) %>% pull()

EN_notHB_asvs <- (setdiff(ENasvs, HBasvs))
EN_notHB_asvs_df <- rfy_plotdat %>% 
  filter(Labeler %in% EN_notHB_asvs) %>% select(ASValias, Reads, order_name, family_name, genus_name, species_name, SiteMonth) %>% 
  group_by(ASValias, order_name, family_name, genus_name, species_name, SiteMonth) %>% 
  summarise(Reads=sum(Reads)) %>% 
  mutate(splitter=ASValias) %>% separate(., splitter, into=c("delete", "ASVnum"), sep = "-") %>% 
  select(-delete) %>% mutate(ASVnum = as.numeric(ASVnum))
HB_notEN_asvs <- (setdiff(HBasvs, ENasvs))
HB_notEN_asvs_df <- rfy_plotdat %>% 
  filter(Labeler %in% HB_notEN_asvs) %>% select(ASValias, Reads, order_name, family_name, genus_name, species_name, SiteMonth) %>% 
  group_by(ASValias, order_name, family_name, genus_name, species_name, SiteMonth) %>% 
  summarise(Reads=sum(Reads)) %>% 
  mutate(splitter=ASValias) %>% separate(., splitter, into=c("delete", "ASVnum"), sep = "-") %>% 
  select(-delete) %>% mutate(ASVnum = as.numeric(ASVnum))


################################################################################
## per ASV summaries
## what's going on at the ASV level... trying to understand why the SL classifier does/doesn't rank an ASV as important 
## helps summarise abundance/occurrence among the ASVs labeled in plot
################################################################################

## 1. Order rank heatmap: all rarefied reads, all ASVs 
ASV_ocr_sumry <- rfy_plotdat %>% 
  group_by(ASValias, SiteMonth) %>% 
  summarise(nSamples = n_distinct(SampleID)) %>% 
  spread(SiteMonth, nSamples, fill = 0) %>% 
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nSamples") %>%
  group_by(SiteMonth) %>%
  mutate(mSamples=sum(nSamples)) %>% 
  mutate(pSamples=round((nSamples/mSamples),2))

ASV_abu_sumry <- rfy_plotdat %>% 
  group_by(ASValias, SiteMonth) %>% 
  summarise(nReads = sum(Reads)) %>% 
  spread(SiteMonth, nReads, fill = 0) %>% 
  gather(`EN-June`, `EN-July`, `EN-September`, `HB-June`, `HB-July`, `HB-September`, key = "SiteMonth", value="nReads") %>%
  group_by(SiteMonth) %>%
  mutate(mReads=sum(nReads)) %>% 
  mutate(pReads=round((nReads/mReads),2))

## group ocurrence and abundance information:
ASV_sumry <- merge(ASV_ocr_sumry, ASV_abu_sumry, by = c('ASValias', 'SiteMonth')) %>% 
  mutate(splitter=SiteMonth) %>% separate(., col=splitter, into=c("Site", "Month"), sep="-")


## which samples were used for Predictions?
## are the ASVs not labeled as important just missing from the training / testing set?
predictions_df <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/SiteMonth_predictions.tsv", delim = "\t", col_names = TRUE)
trainingASVs <- rfy_plotdat %>% 
  filter(!SampleID %in% predictions_df$SampleID) %>% 
  distinct(ASValias)
non_trainingASVs <- rfy_plotdat %>% 
  filter(SampleID %in% predictions_df$SampleID) %>% 
  distinct(ASValias)
## looking at these two data.frames and realizing that the ASVs in our dipteran only plots pretty much are always in both sets
## so they're probably not missing... the classifier training set has most ASVs - so they're just not important because the classifier dosn't think they're important (not because they're missing from the classifier)

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
## big change in analysis below...
## Next set of plots do the same as above for are identifying ASVs trainied for ..
## .. classifier to identify differences in Month WITHOUT ANY SITE class
## Thus observed monthly differences may be due to changes in that Month ..
## .. in one site only, or both sites
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
## create new ASV summary dataset that partitions by Month but aggregates ASV information across both Sites
## pulling select ASVs for "Month" machine learning output
## selecting features in 75th percentile of "Importance" 
filt75_m <- quantile(rmonth_df$importance, .75)
select_MonthASVs <-rmonth_df %>% filter(importance >= filt75_m) %>% pull(ASVid)

## how many of these are the same ASVs? 
length(intersect(select_MonthASVs, selectASVs)) ## 30 are shared (among 60 and 32 for SiteMonth and Month, respectively)
## that's important - the two different models are converging on a similar set of winners

## are any of these actually not listed, or just not in the 75th percentiles?
all_month_importance_ASVs <- rmonth_df$ASVid
length(all_month_importance_ASVs) ## 128 ASVs used in that classifier (about half as many as SiteMonth)
all_sitemonth_importance_ASVs <- rsitemonth_df$ASVid
length(all_sitemonth_importance_ASVs) ## 237 ASVs
intersect(all_month_importance_ASVs, all_sitemonth_importance_ASVs)
## every "Month" ASV is found in "SiteMonth" classifier !!

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

## plot all in a giant plot; save as "slopeplot_allDat_byMonth"; export at 1000x500
ggarrange(ml_sp_ocr_Monthly, ml_sp_abu_Monthly, common.legend = TRUE, labels = c("A", "B"))

## set sample detections and fractions of read abundances as X/Y scatterplot
## save as 'ml_pAbu_by_pOcr_scatterplot_byMonth'; export at 1000x1000

## split up the dataset into three parts:
## 1. ML data (solid circles, colored by 5 Orders now)
## 2. non ML data in same 5 Orders as #1
## 3. non ML data in Orders other than #2

MLasvs_Monthly <- ASV_sumry_SLfiltd_Monthly$ASValias %>% unique(.)
nonML_selectOrders_ASVs_Monthly <- ASV_sumry_Monthly %>%
  filter(order_name %in% OrderKeepers_monthOnly) %>%
  filter(!ASValias %in% ASV_sumry_SLfiltd_Monthly$ASValias) %>%
  pull(ASValias) %>% unique(.)
nonML_otherOrders_ASVs <- ASV_sumry_Monthly %>% filter(!order_name %in% OrderKeepers_monthOnly) %>% pull(ASValias) %>% unique(.)

## use ASV_sumry_Monthly data.frame but rename "order_name" to "plotOrder_name" and subsitute all other Orders not in OrderKeeper to "Other"
ASV_sumry_Monthly$plotOrder_name <- ASV_sumry_Monthly$order_name
xm <- unique(ASV_sumry_Monthly$order_name) ## all Order names
ym <- setdiff(xm, OrderKeepers) ## all "Other" Order names
  ## type all "ym" entries manually in 'pats_m' vector
pats_m <- c('Blattodea|Entomobryomorpha|Hymenoptera|Isopoda|Mecoptera|Mesostigmata|Neuroptera|Odonata|Orthoptera|Poduromorpha|Sarcoptiformes|Strepsiptera|Trombidiformes|Coleoptera|Ephemeroptera|Trichoptera')
## note we've added 3 new Orders to the previous list for SiteMonth (the 3 Orders dropped form the original 8 colored Orders: Coleoptera, Ephemeroptera, and Trichoptera)  
ASV_sumry_Monthly$plotOrder_name <- str_replace_all(ASV_sumry_Monthly$plotOrder_name, pats_m, "Other")  ## this does the substitution

## use ASV_sumry_Monthly data.frame, but add in field that distinguishes ML-labeled ASVs or not w/ TRUE/FALSE statement
ASV_sumry_Monthly <- ASV_sumry_Monthly %>% mutate(MLsample = ASValias %in% MLasvs_Monthly)
## add in filterer for later faceting probs
ASV_sumry_Monthly$MLsample_filter <- ASV_sumry_Monthly$MLsample
ASV_sumry_Monthly$MLsample_filter <- ifelse(ASV_sumry_Monthly$MLsample_filter==TRUE, "SL-feature", "not-SL-feature")

## add in new color to 5-color palette to reflect "Other" now...
pal6 <- c(pal5, "gray60")

## set levels so that "Other" is placed at bottom of legend:
ASV_sumry_Monthly$plotOrder_name <- factor(ASV_sumry_Monthly$plotOrder_name, levels = c(
  "Araneae", "Diptera","Hemiptera","Lepidoptera","Psocodea","Other"))
ASV_sumry_Monthly$Month <- factor(ASV_sumry_Monthly$Month, levels=c("June", "July", "September"))

## plot; save as 'ml_pAbu_by_pOcr_scatterplot_byMonth'; export at 1200x400
## removed labels to avoid artibrary text labeling (we're already discriminating among ML and non-ML sites...)
ggplot() +
  geom_point(data = ASV_sumry_Monthly,
             aes(y = pReads, x=pSamples, color=plotOrder_name, shape=MLsample),
             alpha = 0.8, size=2.5) +
  facet_wrap(~Month) +
  scale_shape_manual(values=c(0,19)) +
  scale_color_manual(values=pal6) +
  scale_x_continuous(expand = expand_scale(add = 0.3)) +
  geom_text_repel(data = ASV_sumry_Monthly %>% filter(pSamples > 0.37),
                  aes(y = pReads, x=pSamples, color=order_name, label=Labeler), 
                  nudge_x = 5, force=5, direction="y", segment.size = 0.2, segment.alpha = 0.5, seed = 42, size=3) +
  labs(y="fraction Reads", x="fraction Samples", color="Arthropod Order", shape="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=12)) +
  guides(color = guide_legend(nrow=1),
         shape=FALSE)

## can make these facets easier to identify differences with animation
## same plot as above, but eliminating Month facet with animation
## ..transition state is Month
mnth.ani <- ggplot() +
  geom_point(data = ASV_sumry_Monthly,
             aes(y = pReads, x=pSamples, color=plotOrder_name, shape=MLsample),
             alpha = 0.8, size=4) +
  scale_shape_manual(values=c(0,19)) +
  scale_color_manual(values=pal6) +
  scale_x_continuous(expand = expand_scale(add = 0.3)) +
  geom_text_repel(data = ASV_sumry_Monthly %>% filter(pSamples > 0.37 | pReads > 0.03),
                  aes(y = pReads, x=pSamples, color=order_name, label=Labeler), 
                  nudge_x = 5, force=5, direction="y", segment.size = 0.2, segment.alpha = 0.5, seed = 42, size=5) +
  labs(y="fraction Reads", x="fraction Samples", color="Arthropod Order", shape="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=12),
        plot.title = element_text(size=20)) +
  guides(color = guide_legend(nrow=2), shape=FALSE) +
  transition_states(Month, transition_length = 1, state_length = 2) + ggtitle('Month: {closest_state}')

## render:
animate(mnth.ani, renderer = gifski_renderer(loop=TRUE), width=600, height = 600)
anim_save("~/Repos/mysosoup/figures/gifs/month_ASVs.gif")

## Plot to focus on Dipterans
## plot; save as ml_pAbu_by_pOcr_scatterplot_byMonth_DipteraOnly; export at 1200x600
ggplot() +
  geom_point(data = ASV_sumry_Monthly %>% filter(order_name=="Diptera"),
             aes(y = pReads, x=pSamples, shape=MLsample),
             alpha = 0.8, size=2.5, color="black") +
  scale_shape_manual(values=c(0,19)) +
  facet_wrap(~ Month, nrow = 3) +
  scale_x_continuous(expand = expand_scale(add = 0.3)) +
  geom_label_repel(data = ASV_sumry_Monthly %>% filter(order_name == "Diptera") %>% filter(MLsample==TRUE) %>% filter(pSamples > 0.3),
                   aes(y = pReads, x=pSamples, label=Labeler),
                   color='black', fill="white", size = 2.25, force = 5, seed = 42, nudge_x = 5, direction = "y",
                   segment.size = 0.2, segment.alpha=0.5, segment.colour = "black") +
  geom_label_repel(data = ASV_sumry_Monthly %>% filter(order_name == "Diptera") %>% filter(MLsample==FALSE) %>% filter(pReads > 0),
                   aes(y = pReads, x=pSamples, label=Labeler),
                   color='black', fill="white", size = 2.25, force = 5, seed = 42, nudge_x = -5, direction = "y",
                   segment.size = 0.2, segment.alpha=0.5, segment.colour = "black") +
  labs(y="fraction Reads", x="fraction Samples", color="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=16),
        plot.title = element_text(size=20, face = "bold")) +
  guides(color = FALSE, shape=FALSE)

## Animation for Dipteran only:
dipAni_monthly <- ggplot() +
  geom_point(data = ASV_sumry_Monthly %>% filter(order_name=="Diptera"),
             aes(y = pReads, x=pSamples, shape=MLsample),
             alpha = 0.8, size=4, color="black") +
  scale_shape_manual(values=c(0,19)) +
  scale_x_continuous(expand = expand_scale(add = 0.33)) +
  geom_label_repel(data = ASV_sumry_Monthly %>% filter(order_name == "Diptera") %>% filter(MLsample_filter=='SL-feature') %>% filter(pSamples > 0.3),
                   aes(y = pReads, x=pSamples, label=Labeler),
                   color='black', fill="white", size = 3, force = 5, seed = 42, nudge_x = 5, direction="y",
                   segment.size = 0.2, segment.alpha=0.5, segment.colour = "black") +
  geom_label_repel(data = ASV_sumry_Monthly %>% filter(order_name == "Diptera") %>% filter(MLsample_filter=='not-SL-feature') %>% filter(pReads > 0),
                   aes(y = pReads, x=pSamples, label=Labeler),
                   color='black', fill=NA, size = 3, force = 5, seed = 42,nudge_x = -5, direction="y",
                   segment.size = 0.3, segment.alpha=0.5, segment.colour = "black") +
  labs(y="fraction Reads", x="fraction Samples", color="") +
  theme_devon() +
  theme(legend.position = "top", legend.text = element_text(size=16), plot.title = element_text(size=22, face = "bold")) +
  guides(color = FALSE, shape=FALSE) +
  transition_states(Month, transition_length = 1, state_length = 2) + ggtitle('Month: {closest_state}')

## render:
animate(dipAni_monthly, renderer = gifski_renderer(loop=TRUE), width=800, height=800)
anim_save("~/Repos/mysosoup/figures/gifs/month_ASVs_DipteraOnly.gif")


################################################################################
## to show that the ASVs that are generating high read counts / observations aren't being ..
## skewed by a single outlier, let's show the distribution of some of the highest ASVs among ..
## dipteran and non-Dipteran ASVs
################################################################################

## get a list of the top counts/abundances for Dipteran taxa only
diptop5 <- alldat %>% filter(order_name=="Diptera") %>% group_by(SampleID) %>% top_n(5, Reads)
diptop5 %>% group_by(ASValias) %>% tally() %>% arrange(-n)
diptop5 %>% group_by(ASValias) %>% summarise(sumReads=sum(Reads), nSamples=n()) %>% arrange(-sumReads)
## 8 of 10 in top counts/reads are overlapping; keeping all 12 and making a list manually from this output
dipASVs = paste("ASV-", c(6,3,10,5,20,4,18,15,17,8,2,11), sep = "")


## repeat to generate a list of the top counts/abundances for all non-Dipteran taxa
nondiptop5 <- alldat %>% filter(order_name!="Diptera") %>% group_by(SampleID) %>% top_n(5, Reads)
nondiptop5 %>% group_by(ASValias) %>% summarise(sumReads=sum(Reads), nSamples=n()) %>% arrange(-sumReads)
nondiptop5 %>% group_by(ASValias) %>% tally() %>% arrange(-n)
## same as above: keeping 12, though not all non-Dipteran Orders represented here
nondipASVs = paste("ASV-", c(1, 13, 9, 14, 19, 25, 7, 58, 12, 48, 26, 30), sep = "")

topASVdf <- alldat %>% filter(ASValias %in% topASVs)
topASVdf$colorlabel <- ifelse(topASVdf$order_name=="Diptera", "Diptera", "nonDiptera")
## set levels for plot
topASVdf$order_name <- factor(topASVdf$order_name, 
                              levels = c("Araneae","Coleoptera","Diptera","Hemiptera","Lepidoptera","Psocodea","Trichoptera"))
topASVdf$ASValias <- factor(topASVdf$ASValias, levels = c(
  'ASV-2', 'ASV-3', 'ASV-4', 'ASV-5', 'ASV-6', 'ASV-8', 'ASV-10', 'ASV-11', 'ASV-15', 'ASV-17', 'ASV-18', 'ASV-20',
  'ASV-1', 'ASV-13', 'ASV-25', 'ASV-7', 'ASV-12', 'ASV-48', 'ASV-9', 'ASV-26', 'ASV-58', 'ASV-14', 'ASV-19', 'ASV-30'))
## keep same color palette scheme from previous styles
pal7 <- c('#3778bf', '#efb435', 'black', '#ff028d', '#825f87', '#d9544d', '#a87900')

## plot; save as 'perSample_topASV_read_and_counts'; export at 
ggplot(topASVdf, 
       aes(x=ASValias, y=Reads, fill=order_name, color=order_name)) +
  geom_point() +
  scale_y_continuous(trans="log2", labels=comma) +
  #geom_violin(alpha=0.5) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_color_manual(values=pal7) +
  scale_fill_manual(values=pal7) +
  facet_wrap(~colorlabel, scales = "free_x") +
  labs(x="", y="log2 sequence counts") +
  theme_devon()

nondiptop5 %>% group_by(ASValias) %>% summarise(sumReads=sum(Reads), nSamples=n()) %>% arrange(-sumReads)
nondiptop5 %>% group_by(ASValias) %>% tally() %>% arrange(-n)




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
