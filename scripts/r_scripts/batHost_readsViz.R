library(tidyverse)
library(reshape2)
library(scales)

theme_devon <- function () { 
  theme_bw(base_size=12, base_family="Avenir") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="transparent", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA)
    )
}

## import bat data summary 
bat_sumry <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/sample_abundancesummaries_wBatHostData.csv")

## reshape dataset for plot
bat_plotdat <- bat_sumry %>% 
  filter(SampleType == "sample") %>% 
  select(SampleID, BatchType, `Myotis lucifugus`, `Myotis sodalis`, `Nycticeius humeralis`, multihit) %>% 
  melt(id.vars = c("SampleID", "BatchType", "multihit"),
       variable.name = "batSpecies", 
       value.name = "Reads")
bat_plotdat$batSpecies <- gsub("Myotis lucifugus", "Myotis\nlucifugus", bat_plotdat$batSpecies)
bat_plotdat$batSpecies <- gsub("Myotis sodalis", "Myotis\nsodalis", bat_plotdat$batSpecies)
bat_plotdat$batSpecies <- gsub("Nycticeius humeralis", "Nycticeius\nhumeralis", bat_plotdat$batSpecies)
bat_plotdat$multihit <- ifelse(bat_plotdat$multihit==TRUE, "multiple", "single")

## set color scale and levels
pal2 <- c("#653700", '#8eab12')
bat_plotdat$BatchType <- factor(bat_plotdat$BatchType, 
                                levels = c("single", "pool"))

## plot; save as 'batHost_ReadsPerSpecies'; export at 500x500
ggplot(data = bat_plotdat %>% filter(!is.na(multihit)),
       aes(x=batSpecies, y=Reads, color=multihit)) +
  geom_jitter(width = 0.3, alpha=0.8) +
  scale_y_continuous(trans="log10", labels=comma) +
  facet_grid (BatchType ~ .) +
  labs(x="", y="host COI sequences", color="Number of Bat Species Per Sample") +
  scale_color_manual(values = pal2) +
  theme_devon() +
  theme(legend.position = "top")


### How many samples contained a particular bat species?
bat_sumry %>% 
  filter(SampleType == "sample") %>% 
  filter(multihit == FALSE) %>% 
  select(SampleID, BatchType, `Myotis lucifugus`, `Myotis sodalis`, `Nycticeius humeralis`, multihit) %>% 
  melt(id.vars = c("SampleID", "BatchType", "multihit"),
       variable.name = "batSpecies", 
       value.name = "Reads") %>% 
  filter(!is.na(Reads)) %>%
  group_by(batSpecies) %>% 
  summarise(counts = n())

## How many samples contained a particular bat species, for only the single samples?
bat_sumry %>% 
  filter(BatchType == "single") %>% 
  filter(SampleType == "sample") %>% 
  filter(multihit == FALSE) %>% 
  select(SampleID, BatchType, `Myotis lucifugus`, `Myotis sodalis`, `Nycticeius humeralis`, multihit) %>% 
  melt(id.vars = c("SampleID", "BatchType", "multihit"),
       variable.name = "batSpecies", 
       value.name = "Reads") %>% 
  filter(!is.na(Reads)) %>%
  group_by(batSpecies) %>%
  summarise(counts = n())
  ## no duplicate data: 137 Indiana Bat, 5 Little Brown, 2 Evening Bat (out of 196 analyzed)

## How many instances are there where a multi-species sample contains Myotis sodalis (and either MYLU or NYHU)?
bat_sumry %>% 
  filter(BatchType == "single") %>% 
  filter(SampleType == "sample") %>% 
  filter(multihit == TRUE) %>% 
  select(SampleID, BatchType, `Myotis lucifugus`, `Myotis sodalis`, `Nycticeius humeralis`, multihit) %>% 
  filter(!is.na(`Myotis sodalis`)) %>% 
  nrow()
  ## 11 instances where there is an Indiana bat sample and something else 

## How many instances where multi-species doesn't include a Myotis sodalis sample?
bat_sumry %>% 
  filter(BatchType == "single") %>% 
  filter(SampleType == "sample") %>% 
  filter(multihit == TRUE) %>% 
  select(SampleID, BatchType, `Myotis lucifugus`, `Myotis sodalis`, `Nycticeius humeralis`, multihit) %>% 
  filter(is.na(`Myotis sodalis`)) %>% 
  nrow()
  ## zero!

### make a list of any single pellet samples that contain ONLY non-sodalis guano
nonMYSO_sampleIDs <- bat_sumry %>% 
  filter(BatchType == "single") %>% 
  filter(SampleType == "sample") %>% 
  filter(multihit == FALSE) %>% 
  select(SampleID, BatchType, `Myotis lucifugus`, `Myotis sodalis`, `Nycticeius humeralis`, multihit) %>% 
  filter(is.na(`Myotis sodalis`)) %>% 
  select(SampleID) %>% pull()
write.table(nonMYSO_sampleIDs, 
          file = "~/github/mysosoup/data/host/nonMYSO_sampleIDlist.txt",
          quote = FALSE, row.names = FALSE, col.names = FALSE)
