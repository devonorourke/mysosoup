#runonce: install.packages("mvabund")
library(tidyverse)
library(mvabund)
library(qiime2R)

## load data:
qza <- read_qza("~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza")

## import metadata 
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv.gz", col_names = TRUE)
meta <- meta %>% 
  select(SampleID, Roost, CollectionMonth, SampleType, Site) %>% 
  rename('Month' = CollectionMonth) %>% 
  filter(SampleType == 'sample')
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$Month <- ifelse(meta$Month == "6", gsub("6", "June", meta$Month), meta$Month)
meta$Month <- ifelse(meta$Month == "7", gsub("7", "July", meta$Month), meta$Month)
meta$Month <- ifelse(meta$Month == "9", gsub("9", "September", meta$Month), meta$Month)
meta$SiteMonth <- paste(meta$Site, meta$Month, sep="-")
meta <- data.frame(meta) %>% filter(SampleID %in% colnames(qza$data))  ## retain only meta data sampleIDs remaining in table:


## doing the same thing with our rarefied data:
## check how abundances per ASV vary (max would be 5200 )
bat_spp <- mvabund(t(qza$data))  ## read in data into their format first!
#boxplot(bat_spp, horizontal=TRUE, las=1, main="Abundance", xlab="", ylab="")  ## only top abundances plotted; and yep, they're variable!
#meanvar.plot(bat_spp) ## wow, do we ever have a linear relationship between abundance variance and mean abundance
bat_mod1 <- manyglm(bat_spp ~ meta$Month*meta$Site, family="negative_binomial")
#summary(bat_mod1)
anova.manyglm()

anova(bat_mod1)
