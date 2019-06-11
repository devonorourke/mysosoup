library(tidyverse)

meta <- read_csv(file="~/Repos/mysosoup/data/mangan_metadata.csv.gz")

colnames(meta)
head(meta[6:15])

qiimemeta <- meta %>%
  select(SampleID, SampleType, BatchType, CollectionMonth, Site, Roost, Source, ContamArea)
colnames(qiimemeta)[1] <- "#Sample ID"

qiimemeta$Site <- ifelse(qiimemeta$Site == "Egner", gsub("Egner", "EN", qiimemeta$Site), qiimemeta$Site)
qiimemeta$Site <- ifelse(qiimemeta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", qiimemeta$Site), qiimemeta$Site)
qiimemeta$CollectionMonth[is.na(qiimemeta$CollectionMonth)] <- "control"
qiimemeta$SiteRoost <- paste(qiimemeta$Site, qiimemeta$Roost, sep="-")
qiimemeta$SiteRoost <- ifelse(qiimemeta$SiteRoost == "control-control", gsub("control-control", "control", qiimemeta$SiteRoost), qiimemeta$SiteRoost)
qiimemeta$CollectionMonth <- ifelse(qiimemeta$CollectionMonth == "6", gsub("6", "June", qiimemeta$CollectionMonth), qiimemeta$CollectionMonth)
qiimemeta$CollectionMonth <- ifelse(qiimemeta$CollectionMonth == "7", gsub("7", "July", qiimemeta$CollectionMonth), qiimemeta$CollectionMonth)
qiimemeta$CollectionMonth <- ifelse(qiimemeta$CollectionMonth == "9", gsub("9", "September", qiimemeta$CollectionMonth), qiimemeta$CollectionMonth)
qiimemeta$ContamArea[is.na(qiimemeta$ContamArea)] <- "isolate"
qiimemeta$SiteMonth <- paste(qiimemeta$Site, qiimemeta$CollectionMonth, sep="-")
qiimemeta$SiteMonth <- ifelse(qiimemeta$SiteMonth == "control-control", gsub("control-control", "control", qiimemeta$SiteMonth), qiimemeta$SiteMonth)
qiimemeta$SiteRoostMonth <- paste(qiimemeta$Site, qiimemeta$Roost, qiimemeta$CollectionMonth, sep="-")
qiimemeta$SiteRoostMonth <- ifelse(qiimemeta$SiteRoostMonth == "control-control-control", gsub("control-control-control", "control", qiimemeta$SiteRoostMonth), qiimemeta$SiteRoostMonth)

write.table(qiimemeta, row.names = FALSE, quote = FALSE, sep = "\t",
            file="~/Repos/mysosoup/data/metadata/qiime_meta.tsv")
