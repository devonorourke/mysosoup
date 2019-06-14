library(tidyverse)
library(qiime2R)
library(vegan)
library(ape)

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

## import metadata 
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/metadata/mangan_metadata.csv", col_names = TRUE)
tinymeta <- meta %>% select(SampleID, Roost, CollectionMonth, SampleType, Site)
tinymeta$Site <- ifelse(tinymeta$Site == "Egner", gsub("Egner", "EN", tinymeta$Site), tinymeta$Site)
tinymeta$Site <- ifelse(tinymeta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", tinymeta$Site), tinymeta$Site)
tinymeta$CollectionMonth[is.na(tinymeta$CollectionMonth)] <- "control"
tinymeta$Labeler <- paste(tinymeta$Site, tinymeta$Roost, sep="-")
tinymeta$Labeler <- ifelse(tinymeta$Labeler == "control-control", gsub("control-control", "control", tinymeta$Labeler), tinymeta$Labeler)
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "6", gsub("6", "June", tinymeta$CollectionMonth), tinymeta$CollectionMonth)
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "7", gsub("7", "July", tinymeta$CollectionMonth), tinymeta$CollectionMonth)
tinymeta$CollectionMonth <- ifelse(tinymeta$CollectionMonth == "9", gsub("9", "September", tinymeta$CollectionMonth), tinymeta$CollectionMonth)


## function import distances and runs Adonis (one per distance measure)
Adonis.function <- function(qzaPath, BetaMethod, FiltType) {
  qza <- read_qza(qzaPath)
  tmp.meta <- data.frame(tinymeta) %>% filter(SampleID %in% labels(qza$data))
  tmp.adonis <- adonis(qza$data ~ CollectionMonth * Site, data = tmp.meta)
  tmp.adonis <- as.data.frame(tmp.adonis$aov.tab)
  tmp.out <- data.frame(tmp.adonis, BetaMethod, FiltType)
  tmp.out$TestGroup <- row.names(tmp.out)
  tmp.out
}

## set file paths for every distance qza file:
noNTCasvs.bc.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.noNTCasvs.bc.beta.qza"
noNTCasvs.ds.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.noNTCasvs.ds.beta.qza"
noNTCasvs.gu.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.noNTCasvs.gu.beta.qza"
noNTCasvs.uu.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.noNTCasvs.uu.beta.qza"
noNTCasvs.wu.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.noNTCasvs.wu.beta.qza"
wNTCasvs.bc.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.wNTCasvs.bc.beta.qza"
wNTCasvs.ds.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.wNTCasvs.ds.beta.qza"
wNTCasvs.gu.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.wNTCasvs.gu.beta.qza"
wNTCasvs.uu.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.wNTCasvs.uu.beta.qza"
wNTCasvs.wu.path="~/Repos/mysosoup/data/qiime_qza/distances/Mangan.wNTCasvs.wu.beta.qza"

## caclulate Adonis values for each:
adonis_noNTC_bc <- Adonis.function(noNTCasvs.bc.path, "bc", "noNTC")
adonis_noNTC_ds <- Adonis.function(noNTCasvs.ds.path, "ds", "noNTC")
adonis_noNTC_gu <- Adonis.function(noNTCasvs.gu.path, "gu", "noNTC")
adonis_noNTC_uu <- Adonis.function(noNTCasvs.uu.path, "uu", "noNTC")
adonis_noNTC_wu <- Adonis.function(noNTCasvs.wu.path, "wu", "noNTC")
adonis_wNTC_bc <- Adonis.function(wNTCasvs.bc.path, "bc", "wNTC")
adonis_wNTC_ds <- Adonis.function(wNTCasvs.ds.path, "ds", "wNTC")
adonis_wNTC_gu <- Adonis.function(wNTCasvs.gu.path, "gu", "wNTC")
adonis_wNTC_uu <- Adonis.function(wNTCasvs.uu.path, "uu", "wNTC")
adonis_wNTC_wu <- Adonis.function(wNTCasvs.wu.path, "wu", "wNTC")

## combine datasets and write as single file:
all_Adonis <- rbind(adonis_noNTC_bc, adonis_noNTC_ds, adonis_noNTC_gu, adonis_noNTC_uu, adonis_noNTC_wu,
                    adonis_wNTC_bc, adonis_wNTC_ds, adonis_wNTC_gu, adonis_wNTC_uu, adonis_wNTC_wu)

write.csv(all_Adonis, file="~/Repos/mysosoup/data/text_tables/adonis_w.or.no_NTC_testing.csv", row.names = FALSE, quote = FALSE)


##############################
