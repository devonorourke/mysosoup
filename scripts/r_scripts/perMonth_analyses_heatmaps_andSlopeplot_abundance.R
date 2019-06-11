library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)

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

## import relevantFeatures data
revfeat <- read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/Learn/optimized/learn-opt-Month/relevantFeatures/importance.tsv",
                      delim = "\t", col_names = TRUE)
colnames(revfeat)[1] <- "ASVid"
revfeat$cumsum <- cumsum(revfeat$importance)
selectASVs <- revfeat %>% filter(cumsum < 0.5) %>% select(ASVid) %>% pull()

## import taxa data, modify ambiguous_taxa as NA
taxa <- read_delim("~/Repos/mysosoup/data/taxonomy/mangan_tax_p97c94.tsv.gz", delim = "\t")
taxa <- taxa %>% separate(., col = Taxon, sep=';',
                          into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>%
  select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
colnames(taxa)[1] <- "ASVid"

## import metadata and add/modify labels for plot
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv.gz", col_names = TRUE) %>%
  select(-RoostDay, -RoostMonth, -RoostYear)
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$CollectionMonth[is.na(meta$CollectionMonth)] <- "control"
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "6", gsub("6", "June", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "7", gsub("7", "July", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "9", gsub("9", "September", meta$CollectionMonth), meta$CollectionMonth)
meta$ContamArea[is.na(meta$ContamArea)] <- "isolate"
meta$SiteMonth <- paste(meta$Site, meta$CollectionMonth, sep="-")
meta$SiteMonth <- ifelse(meta$SiteMonth == "control-control", gsub("control-control", "control", meta$SiteMonth), meta$SiteMonth)
meta$SiteRoostMonth <- paste(meta$Site, meta$Roost, meta$CollectionMonth, sep="-")
meta$SiteRoostMonth <- ifelse(meta$SiteRoostMonth == "control-control-control", gsub("control-control-control", "control", meta$SiteRoostMonth), meta$SiteRoostMonth)
meta <- meta %>% rename(Month=CollectionMonth)

## import read data; filter out samples not included in Machine Learning analysis
## Add ASValias as "ASV-###' where the number reflects per-ASV abundance (highest abundance is first number)
reads <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/mangan.asvtable.long.csv.gz")
tmp1 <- reads %>% group_by(ASVid) %>% summarise(nReads=sum(Reads)) %>% arrange(-nReads) %>%
  mutate(ASValias=paste0("ASV-", row.names(.))) %>% select(-nReads)
reads <- merge(reads, tmp1)
rm(tmp1)

selectSamples <- read_csv('https://github.com/devonorourke/mysosoup/raw/master/data/Learn/Mangan_noBats_famOnly_min9000seqs_noControls_sampleNames.txt',
                          col_names = FALSE) %>% pull()
reads <- reads %>% filter(SampleID %in% selectSamples)

## merege read, and metadata
plotData <- merge(reads, meta) %>% 
  merge(., taxa) %>%
  filter(!is.na(family_name)) %>%
  filter(phylum_name == "Arthropoda") %>%
  select(-PlateNumber, -PlateIndex)

#### ----------- Occurence data, flowing from Order to Family to ASVs
## using 6 colors to label the Order/Family names through to ASV... keeps things consistent when looking at ASV data
## making custom labeler per plot to make this work
## OrderPlotColors are ('windows blue', 'amber', 'gray20', 'swamp', 'dusty purple', 'crimson')
OrderPlotColors <- c('#3778bf', '#feb308', 'gray20', '#698339', '#825f87', '#8c000f')

## generate the Order data
plotSumry_Order_Ocr <- plotData %>% filter(ASVid %in% selectASVs) %>%
  group_by(order_name, Month) %>%
  summarise(nReads = sum(Reads)) %>%
  spread(Month, nReads, fill = 0) %>% 
  gather(June, July, September, key = "Month", value="nReads")


## and plot
plotSumry_Order_Ocr$Month <- factor(plotSumry_Order_Ocr$Month, levels = c("June", "July", "September"))
left <- ggplot(plotSumry_Order_Ocr, aes(x=Month, y=order_name, fill=nReads)) +
  geom_tile(color='gray70') + 
  coord_equal() +
  scale_fill_gradient(low = 'gray90', high = '#7f5e00', breaks = c(0, 500000), limits =c(0, 500000), labels = comma) +
  labs(x="", y="", fill="Reads") +
  scale_y_discrete(position = "right") +
  theme_devon() +
  theme(legend.position = "top",
        axis.text.y = element_text(color=OrderPlotColors),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## generate the Family data
plotSumry_Family_Ocr <- plotData %>% filter(ASVid %in% selectASVs) %>%
  group_by(Month, order_name, family_name) %>%
  summarise(nReads = sum(Reads)) %>%
  spread(Month, nReads, fill = 0) %>% 
  gather(June, July, September, key = "Month", value="nReads") %>%
  mutate(Labeler = paste(order_name, family_name, sep = " - "))

## need to order Family names by Order first
Family_plotNames <- unique(plotSumry_Family_Ocr$Labeler) %>% gsub("^.* - ", "", .)
plotSumry_Family_Ocr$Month <- factor(plotSumry_Family_Ocr$Month, levels = c("June", "July", "September"))
## custom colors:
FamilyPlotColors <- c('#3778bf', '#feb308', rep('gray20', 5), rep('#698339', 2), rep('#825f87', 2), '#8c000f')

## and plot; export at width=600
center <- ggplot(plotSumry_Family_Ocr, aes(x=Month, y=factor(family_name, level = Family_plotNames), fill=nReads)) +
  geom_tile(color='gray70') + 
  coord_equal() +
  scale_fill_gradient(low = 'gray90', high = '#7f5e00', breaks = c(0, 500000), limits =c(0, 500000), labels = comma) +
  scale_y_discrete(position = "right") +
  labs(x="", y="", fill="Samples")+
  theme_devon() +
  theme(legend.position = "none",
        axis.text.y = element_text(color = FamilyPlotColors),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## generate the ASV data
plotSumry_ASV_Ocr <- plotData %>% filter(ASVid %in% selectASVs) %>%
  group_by(Month, order_name, family_name, ASVid, ASValias) %>%
  summarise(nReads = sum(Reads)) %>%
  spread(Month, nReads, fill = 0) %>% 
  gather(June, July, September, key = "Month", value="nReads")
plotSumry_ASV_Ocr <- plotSumry_ASV_Ocr %>% mutate(Labeler = paste(order_name, family_name, ASValias, sep = " - "))
plotSumry_ASV_Ocr <- merge(plotSumry_ASV_Ocr, taxa) ## adding and modifying genus/species info for ASV plot
plotSumry_ASV_Ocr$species_name <- plotSumry_ASV_Ocr$species_name %>% as.character(.) %>% replace_na("unassigned")
plotSumry_ASV_Ocr$genus_name <- plotSumry_ASV_Ocr$genus_name %>% as.character(.) %>% replace_na("unassigned")
plotSumry_ASV_Ocr$labelName <- paste(plotSumry_ASV_Ocr$genus_name, plotSumry_ASV_Ocr$species_name, sep = " ") ## use this label for plot

## order by Order and Family 
ASV_plotNames <- sort(plotSumry_ASV_Ocr$Labeler) %>% gsub("^.* - ", "", .)

## plot - note we're using a slopegraph here instead of a heatmap...
#!orderColors <- c('crimson', 'puce', 'light grey', 'moss green', 'sea blue', 'goldenrod')
SelectorderColors <- c('#3778bf', '#feb308', '#698339', '#825f87', '#8c000f')

## generate left/right labels for plot
leftASVs <- plotSumry_ASV_Ocr %>% filter(order_name=="Diptera") %>% distinct(ASVid) %>% pull(ASVid)
rightASVs <- plotSumry_ASV_Ocr %>% filter(order_name!="Diptera") %>% distinct(ASVid) %>% pull(ASVid)

plotSumry_ASV_Ocr$Month <- factor(plotSumry_ASV_Ocr$Month, levels = c("June", "July", "September"))
right <- ggplot(data = plotSumry_ASV_Ocr,
                aes(x = Month, y = nReads, group = ASVid, label = ASValias)) +
  scale_color_manual(values=orderColors) +
  scale_y_continuous(labels = comma) +
  geom_line(data = plotSumry_ASV_Ocr %>% filter(order_name == "Diptera"), color="grey75", size = 1.25) +
  geom_line(data = plotSumry_ASV_Ocr %>% filter(order_name != "Diptera"), aes(color=order_name), size = 1.25) +
  geom_label_repel(data = plotSumry_ASV_Ocr %>% filter(Month=="June" & ASVid %in% leftASVs & nReads >= 1), 
                   color = "grey30", fill="white",
                   nudge_x = -5, direction = "y", size=3, segment.size = 0.2, segment.colour = "grey70") +
  geom_label_repel(data = plotSumry_ASV_Ocr %>% filter(Month=="September" & ASVid %in% rightASVs), 
                   aes(color = order_name),
                   nudge_x = 5, direction = "y", size=3, segment.size = 0.2, segment.colour = "grey70") +
  labs(x="", y = "Read Abundance", color = "Order") +
  theme_devon() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


plot_grid(left, center, right, ncol = 3, rel_widths = c(1, 1, 3))
