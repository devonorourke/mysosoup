## script creates heatmaps of Machine learning accuracy output from QIIME2 'classify samples'

library(tidyverse)
library(reshape2)

## function for plot theme:
theme_devon <- function () {theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(panel.background  = element_blank(), plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA), legend.key = element_rect(fill="transparent", colour=NA))}


## going to repeat this plot for both SiteMonth and for Month

################################################################################
## plot 1 - Site Month
################################################################################

## import data (accuracy data and prediction sample sizes are separate files):
mlSiteMonth_raw <- data.frame(read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/SiteMonth_accuracy.tsv", delim = "\t", col_names = TRUE))
row.names(mlSiteMonth_raw) <- mlSiteMonth_raw$X1
mlSiteMonth_df <- mlSiteMonth_raw[c(1:6), c(1:7)] %>% 
  melt(., id.vars = c("X1")) %>% 
  rename(X1 = 'actual', variable = 'predicted')

mlSiteMonth_predict <- data.frame(read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/SiteMonth_predictions.tsv", delim= "\t", col_names=TRUE)) %>% 
  group_by(prediction) %>% 
  summarise(counts=n()) %>% 
  rename(prediction='SiteMonth')

mlSiteMonth_plot <- merge(mlSiteMonth_df, mlSiteMonth_predict, by.x='actual', by.y='SiteMonth') %>% 
  mutate(xLabeler = paste(actual, paste("(n = ", counts, ")", sep=""), sep = "\n"))

## sub the '.' in the y-axis labels to match xaxis '-' for labels
mlSiteMonth_plot$predicted <- gsub("\\.", "\n", mlSiteMonth_plot$predicted)

## set levels for plot
mlSiteMonth_plot$predicted <- factor(mlSiteMonth_plot$predicted, levels = c("EN\nJune", "EN\nJuly", "EN\nSeptember", "HB\nJune", "HB\nJuly", "HB\nSeptember"))
mlSiteMonth_plot$xLabeler <- factor(mlSiteMonth_plot$xLabeler, levels = c(
  "EN-June\n(n = 7)", "EN-July\n(n = 11)", "EN-September\n(n = 7)", 
  "HB-June\n(n = 11)", "HB-July\n(n = 9)", "HB-September\n(n = 11)"))

## and plot; save as 'ml_heatmap_SiteMonth', export at 550x450
ggplot(mlSiteMonth_plot, aes(x=xLabeler, y=predicted, fill=value)) +
  geom_tile(stat="identity") +
  coord_flip() +
  scale_fill_gradient(low="#fff8dc", high="#512698") +
  labs(x="true group", y="predicted group", 
       fill="proportion", caption = paste(paste("Overall accuracy: ", round(mlSiteMonth_raw[7,8],2), sep=" "),
      paste("Accuracy ratio: ", round(mlSiteMonth_raw[9,8],2), sep = " "),
       sep = '\n')) +
  theme_devon()


################################################################################
## plot 2 - Month only
################################################################################

## import accuracy data:
mlMonth_raw <- data.frame(read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/Month_accuracy.tsv", delim = "\t", col_names = TRUE))
row.names(mlMonth_raw) <- mlMonth_raw$X1
mlMonth_df <- mlMonth_raw[c(1:3), c(1:4)] %>% 
  melt(., id.vars = c("X1")) %>% 
  rename(X1 = 'actual', variable = 'predicted')

mlMonth_predict <- data.frame(read_delim(file="https://github.com/devonorourke/mysosoup/raw/master/data/MachineLearn/rarefy/Month_predictions.tsv", delim= "\t", col_names=TRUE)) %>% 
  group_by(prediction) %>% 
  summarise(counts=n()) %>% 
  rename(prediction='group')

mlMonth_plot <- merge(mlMonth_df, mlMonth_predict, by.x='actual', by.y='group') %>% 
  mutate(xLabeler = paste(actual, paste("(n = ", counts, ")", sep=""), sep = "\n"))


## set levels for plot
mlMonth_plot$predicted <- factor(mlMonth_plot$predicted, levels = c("June", "July", "September"))
mlMonth_plot$xLabeler <- factor(mlMonth_plot$xLabeler, levels = c(
  "June\n(n = 22)", "July\n(n = 16)", "September\n(n = 18)"))

## and plot; save as 'ml_heatmap_Month', export at 550x425
ggplot(mlMonth_plot, aes(x=xLabeler, y=predicted, fill=value)) +
  geom_tile(stat="identity") +
  coord_flip() +
  scale_fill_gradient(low="#fff8dc", high="#512698") +
  labs(x="true group", y="predicted group", 
       fill="proportion", caption = paste(paste("Overall accuracy: ", round(mlMonth_raw[4,5],2), sep=" "),
                                          paste("Accuracy ratio: ", round(mlMonth_raw[6,5],2), sep = " "),
                                          sep = '\n')) +
  theme_devon()
