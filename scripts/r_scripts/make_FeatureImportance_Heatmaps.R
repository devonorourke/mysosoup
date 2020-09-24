## script creates heatmaps of Machine learning accuracy output from QIIME2 'classify samples'

library(tidyverse)
library(scico)
library(ggpubr)
library(qiime2R)


################################################################################
## step 1 - import data and reformat for plotting
################################################################################

## import taxonomy information
taxa <- read_csv(file="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/taxonomy/filtd_tax_dataframe_ALL.csv")
colnames(taxa)[1] <- "OTUid"

## function to create output for plots
reformaterFunc <- function(featPath, Labeler){
  imp_tmp <- data.frame(read_delim(file=featPath, delim = "\t", skip = 1, col_names = FALSE)) %>% 
    rename(OTUid = X1, Value = X2) %>% 
    mutate(CumValue = cumsum(Value))
  merge(imp_tmp, taxa) %>% 
    mutate(GroupLabel = Labeler) %>% 
    arrange(CumValue)
}

## paths to data to import
Month_featPath="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/MachineLearn/month_importance.tsv"
Site_featPath="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/MachineLearn/site_importance.tsv"
SiteMonth_featPath="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/MachineLearn/sitemonth_importance.tsv"

## import data with 'reformaterFunc' function
Month_featImp_df <- reformaterFunc(Month_featPath, "Month")
Site_featImp_df <- reformaterFunc(Site_featPath, "Site")
SiteMonth_featImp_df <- reformaterFunc(SiteMonth_featPath, "SiteMonth")

rm(taxa)

## combine datasets:
all_featImp_df <- rbind(Month_featImp_df, Site_featImp_df, SiteMonth_featImp_df)

################################################################################
## step 2a - plot each OTU relative importance as heatmap.
################################################################################

## generate color palette to use for arthropod orders:
#original design modified from:  scico(7, palette = 'batlow')
#order8pal <- c("#001959", "#F9CCF9", "tan4", "#3E6C54", "#FDAC9F", "#808133", "#D49347", "turquoise4")
## becuase we'll order the colors alphabetically (by arthropod order), these should be:
## "#001959"  == Araneae
## "#F9CCF9" == Coleoptera
## "#tan4"  == Diptera
## "#3E6C54" == Ephemeroptera
## "#FDAC9F" == Hemiptera
## "#808133" == Lepidoptera
## "#D49347" == Psocodoea
## "turquoise4" == Trichoptera

## reorder the plot by Decending abundance values
## mutate the labels to include family/genus names if missing genus/species names
all_featImp_df <- all_featImp_df %>% 
  mutate(OTUid = reorder(OTUid, -Value)) %>% 
  mutate(Genus = ifelse(is.na(Genus), paste0("Family: ",Family), Genus)) %>% 
  mutate(Species=ifelse(is.na(Species), paste0("Genus: ", Genus), Species)) %>% 
  mutate(Species=gsub("Genus: Family:", "Family:", Species))
## add in a label name for each X axis OTU to include the special name, along with the generic OTU#:
PlotLabel_df <- all_featImp_df %>% 
  select(OTUid, Species) %>% 
  distinct() %>% 
  mutate(PlotName = paste0("OTU-", rownames(.), ": ", Species))
all_featImp_df <- merge(all_featImp_df, PlotLabel_df, by=c('OTUid', 'Species'))

## regroup the PlotName positions by arthropod Order, maintaining initial desecending value ordering
all_featImp_df$PlotName <- fct_reorder(all_featImp_df$PlotName, all_featImp_df$Order, max)

# ## plot
# ggplot(all_featImp_df, 
#        aes(x = PlotName, 
#            y = GroupLabel, 
#            fill=Value)) +
#   geom_tile() +
#   coord_equal() +
#   scale_fill_scico(palette = "lajolla", limits = c(0, 0.15), breaks=c(0, 0.05, 0.10, 0.15)) +
#   labs(x="", y="", fill="Relative importance") +
#   theme(panel.background = element_blank(),
#         legend.position="top",
#         axis.text.x = element_text(angle=90, hjust=1, size=8))
# 
# ## export this plot, but modify x axis labels to color by 8 palette scheme proposed above
# ggsave("FigureX_RelativeFeatureImportance_CoreFeatures.png", height=10, width = 35, units="cm")
# ggsave("FigureX_RelativeFeatureImportance_CoreFeatures.pdf", height=10, width = 35, units="cm")


## could also split up plots by faceting using taxonomic order:
## plot
p1 <- ggplot(all_featImp_df, 
       aes(x = PlotName, 
           y = GroupLabel, 
           fill=Value)) +
  geom_tile() +
  #coord_equal() +
  facet_grid(~Order, scales = "free_x", space = "free_x") +
  #scale_fill_scico(palette = "broc", direction = 1, begin = 0.49,
  scale_fill_scico(palette = "cork", direction = 1, begin = 0.51,
                   limits = c(0, 0.15), breaks=c(0, 0.05, 0.10, 0.15)) +
  labs(x="", y="", fill="Relative importance") +
  theme(legend.position="bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, size=8),
        strip.text = element_blank(),
        panel.spacing = unit(1, "lines"))

################################################################################
## alternative plot to add in abundance information?
################################################################################

## add in read abundance info:

## download from github:
#download.file(url = 'https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/Mangan.clust_p985_table_Filtd_min10k.qza')
#coretable <- read_qza('Mangan.clust_p985_table_Filtd_min10k.qza')

## or run from local:
coretable <- read_qza('~/github/mysosoup/data/qiime_qza/Mangan.clust_p985_table_Filtd_min10k.qza')
coredata <- as.data.frame(coretable$data) %>% 
  mutate(OTUid = row.names(.)) %>% 
  pivot_longer(-OTUid, names_to = "SampleID", values_to = "Reads") %>% 
  filter(Reads > 0)
rm(coretable)

coresamplesumry <- coredata %>% 
  group_by(SampleID) %>% 
  summarise(TotalReads = sum(Reads))

readdata <- merge(coredata, coresamplesumry, by='SampleID') %>% 
  mutate(FracReads = Reads / TotalReads) %>% 
  mutate(LogFracReads = log10(FracReads)) %>% 
  mutate(LogFracReads = round(LogFracReads, 3))

extrataxadata <- all_featImp_df %>% select(OTUid, Order, PlotName) %>% distinct()
readdata <- merge(readdata, extrataxadata, by="OTUid")

meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", 
                 col_names = TRUE) %>% 
  filter(SampleType == "sample") %>% 
  select(SampleID, CollectionMonth, Site) %>% 
  mutate(Site = case_when(Site == "Egner" ~ "EN", Site == "HickoryBottoms" ~ "HB")) %>% 
  mutate(Month = case_when(CollectionMonth == 6 ~ "June", CollectionMonth == 7 ~ "July", CollectionMonth == 9 ~ "Sept")) %>% 
  select(-CollectionMonth)

readdata <- merge(readdata, meta, all.x=TRUE, by="SampleID")

## regroup the PlotName positions by arthropod Order, maintaining initial desecending value ordering
readdata$PlotName <- fct_reorder(readdata$PlotName, readdata$Order, max)

## set color palette:
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)

## order the Months:
readdata$Month <- factor(readdata$Month, levels = c("June", "July", "Sept"))

p2 <- ggplot(readdata, 
       aes(x=PlotName,
           y=FracReads,
           shape=Site,
           color=Month)) +
  geom_point() +
  scale_color_manual(values=v3pal) +
  facet_grid(~Order, scales = "free_x", space = "free_x") +
  labs(x="", y="Fraction reads per sample") +
  theme(legend.position = "top",
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        #strip.text = element_blank(),
        panel.spacing = unit(1, "lines"))

## paste together:
ggarrange(p2, p1, nrow = 2, align="v", heights = c(0.6, 1), labels=c("D"))

ggsave("~/github/mysosoup/figures/FigureX_RelativeFeatureImportance_CoreFeatures_wAbundances.png", height=20, width = 35, units="cm")
ggsave("~/github/mysosoup/figures/FigureX_RelativeFeatureImportance_CoreFeatures_wAbundances.pdf", height=29, width = 35, units="cm")

## modify .pdf to include silhouettes for publication

###########################3

## select a subset of these OTUs to get further info about taxa:
selectOTUs <- c('65d5fa6d7e70a2048699fd898caf1fca',
                '954ba86f46a17e332fb8f7e065d5c196',
                'f1c7baf332ebaa98a4381b5c5de74f64',
                '5cac244b7bc2b762706cf1e5898f59a9',
                '2066a2ef9e7c8048b1af441a0d9614b1',
                'ee428b867f1c0aa2d56e7c61d744bf56',
                '6a5c7f2f53495aec6dc9eea864305059',
                '46f68fb7accd2cf2419ca3cd22da3646',
                '2c11ea64de16805e25ed359888baabf6',
                '6280a664e750c88afb6069f05b2a29a7',
                '9c1759410185e12aa041d91bef2484d5')

## subset the data from these select ones... this is kinda arbritrary 'eyeballing' ones though...
selectOTUdat <- readdata %>% 
  filter(OTUid %in% selectOTUs) %>% 
  group_by(OTUid, PlotName, Site, Month) %>% 
  summarise(nReads=sum(Reads), nDetections=n())

## what about filtering using the most Importance to each group?
## get the cummulative sums for each group (site, date, site-date), and sort out:
tmp <- all_featImp_df %>% select(-Classifier, -Confidence, -Class, -Consensus) %>% arrange(Value)
tmp2 <- tmp %>% filter(Value > 0.005) %>% group_by(GroupLabel) %>% select(PlotName, Value, GroupLabel) %>% pivot_wider(names_from = "GroupLabel", values_from = "Value")
tmp3 <- tmp2 %>% select(PlotName, SiteMonth) %>% filter(!is.na(SiteMonth)) %>% arrange(-SiteMonth) %>% mutate(CumImpSiteMonth = cumsum(SiteMonth)) %>% select(-SiteMonth)
tmp4 <- tmp2 %>% select(PlotName, Site) %>% filter(!is.na(Site)) %>% arrange(-Site) %>% mutate(CumImpSite = cumsum(Site)) %>% select(-Site)
tmp5 <- tmp2 %>% select(PlotName, Month) %>% filter(!is.na(Month)) %>% arrange(-Month) %>% mutate(CumImpMonth = cumsum(Month)) %>% select(-Month)
tmp6 <- merge(tmp3, tmp4, by = "PlotName", all=TRUE)
cumImpdf <- merge(tmp6, tmp5, by = "PlotName", all=TRUE)
write_csv(cumImpdf, path="~/github/mysosoup/data/MachineLearn/cumImportance_allGroups.csv")

## now get the taxa for the OTUs that provided the top 50% importance among these 3 groups:
top50Imp_OTUnames <- cumImpdf %>% 
  pivot_longer(-PlotName, names_to = "Group") %>% 
  filter(value <= 0.501) %>% 
  pivot_wider(names_from = "Group", values_from = "value") %>% 
  select(PlotName) %>% pull()

## and reorgainze as before, but with this select list:
top50Imp_dat <- readdata %>% 
  filter(PlotName %in% top50Imp_OTUnames) %>% 
  group_by(OTUid, PlotName, Site, Month) %>% 
  summarise(nReads=sum(Reads), nDetections=n())

