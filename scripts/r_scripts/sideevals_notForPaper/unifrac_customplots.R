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


library(qiime2R)
library(tidyverse)
library(vegan)
library(ape)

## import distances from QIIME
unwt_qza <- read_qza(file="~/Repos/mysosoup/data/qiime_qza/distances/mangan_nobats_uuni_dist.qza")
unwt_dist <- unwt_qza$data

wuf_qza <- read_qza(file="~/Repos/mysosoup/data/qiime_qza/distances/mangan_nobats_wuni_dist.qza")
wuf_dist <- wuf_qza$data

## import select metadata:
meta <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/mangan_metadata.csv.gz", col_names = TRUE)
tinymeta <- meta %>% select(SampleID, Roost, CollectionMonth, Site)
tinymeta$Site <- ifelse(tinymeta$Site == "Egner", gsub("Egner", "EN", tinymeta$Site), tinymeta$Site)
tinymeta$Site <- ifelse(tinymeta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", tinymeta$Site), tinymeta$Site)
tinymeta$CollectionMonth[is.na(tinymeta$CollectionMonth)] <- "control"
tinymeta$Labeler <- paste(tinymeta$Site, tinymeta$Roost, sep="-")
tinymeta$Labeler <- ifelse(tinymeta$Labeler == "control-control", gsub("control-control", "control", tinymeta$Labeler), tinymeta$Labeler)
tinymeta$CollectionMonth <- as.factor(tinymeta$CollectionMonth)
rm(meta)

## run PCoA
pcoa.uuf <- pcoa(unwt_dist)
pcoa.wuf <- pcoa(wuf_dist)

## gather data to plot
pcoa.uuf.vec <- data.frame(pcoa.uuf$vectors) %>% mutate(SampleID = row.names(.))
pcoa.wuf.vec <- data.frame(pcoa.wuf$vectors) %>% mutate(SampleID = row.names(.))

## merge with metadata
pcoa.uuf_df <- merge(pcoa.uuf.vec, tinymeta)
pcoa.wuf_df <- merge(pcoa.wuf.vec, tinymeta)
rm(pcoa.uuf.vec, pcoa.wuf.vec)

## gather fraction of variance explained per Axis (for plot labels)
pcoa.uuf.eig <- data.frame(pcoa.uuf$values) %>% 
  select(Cum_corr_eig) %>% 
  mutate(EigVec = paste0("Axis.", row.names(.)))
pcoa.uuf.eig <- pcoa.uuf.eig %>% 
  mutate(FracVar = c(pcoa.uuf.eig$Cum_corr_eig[1],diff(pcoa.uuf.eig$Cum_corr_eig))) %>%
  select(EigVec, FracVar)

pcoa.wuf.eig <- data.frame(pcoa.wuf$values) %>% 
  select(Cum_corr_eig) %>% 
  mutate(EigVec = paste0("Axis.", row.names(.))) 
pcoa.wuf.eig <- pcoa.wuf.eig %>%
  mutate(FracVar = c(pcoa.wuf.eig$Cum_corr_eig[1],diff(pcoa.wuf.eig$Cum_corr_eig))) %>%
  select(EigVec, FracVar)

## plot
uuf.xAxis_val <- pcoa.uuf.eig %>% filter(EigVec == "Axis.1") %>% select(FracVar) %>% pull()
uuf.xAxis_val <- round(uuf.xAxis_val*100, 2)
uuf.yAxis_val <- pcoa.uuf.eig %>% filter(EigVec == "Axis.2") %>% select(FracVar) %>% pull()
uuf.yAxis_val <- round(uuf.yAxis_val*100, 2)

wuf.xAxis_val <- pcoa.wuf.eig %>% filter(EigVec == "Axis.1") %>% select(FracVar) %>% pull()
wuf.xAxis_val <- round(wuf.xAxis_val*100, 2)
wuf.yAxis_val <- pcoa.wuf.eig %>% filter(EigVec == "Axis.2") %>% select(FracVar) %>% pull()
wuf.yAxis_val <- round(wuf.yAxis_val*100, 2)

## unweighted:
ggplot(pcoa.uuf_df, aes(x=Axis.1, y=Axis.2, color=CollectionMonth, shape=Site)) + 
  geom_point() +
  theme_devon() + 
  labs(x=paste("Axis.1  [", uuf.xAxis_val, "%]", sep = ""),
       y=paste("Axis.2  [", uuf.yAxis_val, "%]", sep = ""),
       caption = "Unweighted Unifrac distances")

## weighted
ggplot(pcoa.wuf_df, aes(x=Axis.1, y=Axis.2, color=CollectionMonth, shape=Site)) + 
  geom_point() +
  theme_devon() + 
  labs(x=paste("Axis.1  [", wuf.xAxis_val, "%]", sep = ""),
       y=paste("Axis.2  [", wuf.yAxis_val, "%]", sep = ""),
       caption = "Unweighted Unifrac distances")

## ---------------
pca.uuf <- data.frame(cmdscale(unwt_dist, k=2)) %>% mutate(SampleID=row.names(.))
ggplot(pca.uuf, aes(x=X1, y=X2)) + geom_point()
