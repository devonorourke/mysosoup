## plot for Figure 4 of manuscript
## showing accuracy of model for Site, Month, and SiteMonth classifiers

library(tidyverse)
library(scico)
library(ggpubr)
library(svglite)

################################################################################
## part 1 - import data and reformat for plots
################################################################################

meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", 
                 col_names = TRUE) %>% 
  filter(SampleType == "sample") %>% 
  select(SampleID, CollectionMonth, Site) %>% 
  mutate(Site = case_when(Site == "Egner" ~ "EN", Site == "HickoryBottoms" ~ "HB")) %>% 
  mutate(Month = case_when(CollectionMonth == 6 ~ "June", CollectionMonth == 7 ~ "July", CollectionMonth == 9 ~ "Sept")) %>% 
  mutate(SiteMonth = paste0(Site, "-", Month)) %>% 
  select(-CollectionMonth)

## for CollectionMonth data:
MonthURL="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/MachineLearn/month_predictions.tsv"
month_pred <- read_delim(file=MonthURL, delim="\t", col_names = TRUE) %>% rename(Prediction = prediction) %>% mutate(Prediction = ifelse(Prediction == "September", "Sept", Prediction))
month_dat <- merge(month_pred, meta, by="SampleID", all.x = TRUE)
month_counts <- month_dat %>% group_by(Month) %>% summarise(ActualCounts = n())
month_sums <- month_dat %>% group_by(Prediction, Month) %>% summarise(Counts = n())
month_plotdat <- merge(month_sums, month_counts) %>% mutate(FractionCounts = Counts/ActualCounts)
rm(month_pred, month_dat, month_counts, month_sums, MonthURL)

## for Site data:
SiteURL="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/MachineLearn/site_predictions.tsv"
site_pred <- read_delim(file=SiteURL, delim="\t", col_names = TRUE) %>% rename(Prediction = prediction)
site_dat <- merge(site_pred, meta, by="SampleID", all.x = TRUE)
site_counts <- site_dat %>% group_by(Site) %>% summarise(ActualCounts = n())
site_sums <- site_dat %>% group_by(Prediction, Site) %>% summarise(Counts = n())
site_plotdat <- merge(site_sums, site_counts) %>% mutate(FractionCounts = Counts/ActualCounts)
rm(site_pred, site_dat, site_counts, site_sums, SiteURL)

## for SiteMonth data
SiteMonthURL = "https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/MachineLearn/sitemonth_predictions.tsv"
sitemonth_pred <- read_delim(file=SiteMonthURL, delim="\t", col_names = TRUE) %>% rename(Prediction = prediction) %>% 
  mutate(Prediction = gsub("September", "Sept", Prediction))
sitemonth_dat <- merge(sitemonth_pred, meta, by="SampleID", all.x = TRUE)
sitemonth_counts <- sitemonth_dat %>% group_by(SiteMonth) %>% summarise(ActualCounts = n())
sitemonth_sums <- sitemonth_dat %>% group_by(Prediction, SiteMonth) %>% summarise(Counts = n())
sitemonth_plotdat <- merge(sitemonth_sums, sitemonth_counts) %>% mutate(FractionCounts = Counts/ActualCounts)
rm(sitemonth_pred, sitemonth_dat, sitemonth_counts, sitemonth_sums, SiteMonthURL)

rm(meta)

################################################################################
## part 2 - make individual plots
################################################################################

## Month plot
## set levels for plot:
month_plotdat$Month <- factor(month_plotdat$Month, levels = c("June", "July", "Sept"))
month_plotdat$Prediction <- factor(month_plotdat$Prediction, levels = c("June", "July", "Sept"))
## and plot
p1 <- ggplot(month_plotdat, aes(x=Prediction, y=Month, fill=FractionCounts, label=Counts)) +
  geom_tile() +
  geom_text(data = month_plotdat %>% filter(Counts < 20), size=8) +
  geom_text(data = month_plotdat %>% filter(Counts >= 20), size=8, color="white") +
  scale_fill_scico(palette = "bilbao", direction = 1, breaks = c(0.25, 0.5, 0.75)) +
  theme_bw() +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        panel.grid = element_blank()) +
  labs(fill = "Proportion of samples\naccurately classified",
       x="\nPredicted", y="Actual\n")

p1
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_byMonth.png", height = 9.5, width = 12, units="cm")
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_byMonth.pdf", height = 9.5, width = 12, units="cm")
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_byMonth.svg", height = 9.5, width = 12, units="cm")

## Site plot
p2 <- ggplot(site_plotdat, aes(x=Prediction, y=Site, fill=FractionCounts, label=Counts)) +
  geom_tile() +
  geom_text(data = site_plotdat %>% filter(Counts < 30), size=8) +
  geom_text(data = site_plotdat %>% filter(Counts >= 30), size=8, color="white") +
  scale_fill_scico(palette = "bilbao", begin=0, end=1, direction = 1, breaks = c(0.25, 0.5, 0.75)) +
  theme_bw() +
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=14),
        panel.grid = element_blank()) +
  labs(fill = "Proportion of samples\naccurately classified",
       x="\nPredicted", y="Actual\n")

p2
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_bySite.png", height = 9.5, width = 12, units="cm")
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_bySite.pdf", height = 9.5, width = 12, units="cm")
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_bySite.svg", height = 9.5, width = 12, units="cm")

## SiteMonth plot
## set levels for plot
sitemonth_plotdat$SiteMonth <- factor(sitemonth_plotdat$SiteMonth, levels=c(
  "EN-June", "EN-July", "EN-Sept", "HB-June", "HB-July", "HB-Sept"))
sitemonth_plotdat$Prediction <- factor(sitemonth_plotdat$Prediction, levels=c(
  "EN-June", "EN-July", "EN-Sept", "HB-June", "HB-July", "HB-Sept"))
## and plot
p3 <- ggplot(sitemonth_plotdat, aes(x=Prediction, y=SiteMonth, fill=FractionCounts, label=Counts)) +
  geom_tile() +
  geom_text(data = sitemonth_plotdat %>% filter(Counts < 20), size=8) +
  geom_text(data = sitemonth_plotdat %>% filter(Counts >= 20), size=8, color="white") +
  scale_fill_scico(palette = "bilbao", direction = 1, breaks = c(0.25, 0.5, 0.75)) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title = element_text(size=14),
        panel.grid = element_blank()) +
  labs(fill = "Proportion of samples\naccurately classified",
       x="\nPredicted", y="Actual\n")

p3
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_bySiteMonth.png", height = 9.5, width = 12, units="cm")
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_bySiteMonth.pdf", height = 9.5, width = 12, units="cm")
ggsave("~/github/mysosoup/figures/FigureY_MLheatmap_bySiteMonth.svg", height = 9.5, width = 12, units="cm")

##Stitch all 3 plots togehter?
p1a <- p1 + labs(x="", y="Actual")
p2a <- p2 + labs(x="Predicted", y="")
p3a <- p3 + labs(x="", y="")
ggarrange(p1a, p2a, p3a, common.legend = TRUE, ncol=3, align = "h", labels = c("A", "B", "C"))

ggsave("~/github/mysosoup/figures/FigureY_MLheatmaps_all.png", height=12, width = 30, units="cm")
ggsave("~/github/mysosoup/figures/FigureY_MLheatmaps_all.pdf", height=12, width = 30, units="cm")
ggsave("~/github/mysosoup/figures/FigureY_MLheatmaps_all.svg", height=12, width = 30, units="cm")
