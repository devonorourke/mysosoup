#runonce: library(devtools)
#runonce: devtools::install_github("benjjneb/decontam")
library(phyloseq)
library(decontam)
library(tidyverse)
library(scales)

sumry0 <- as.data.frame(df %>% group_by(SampleID, SampleType, BatchType) %>% 
  summarise(sumReads=sum(Reads), nASVs=n_distinct(ASVid)) %>%
  arrange(sumReads))

sumry0$index <- seq(nrow(sumry0))
sumry0$is.neg <- sumry0$SampleType == "control"
row.names(sumry0) <- sumry0$SampleID

## plot number of reads per sample; color by $SampleType
ggplot(sumry0, aes(x=index, y=sumReads, color=SampleType)) + 
  geom_point() +
  geom_point(data=sumry0 %>% filter(SampleType=="control"), size=3) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = comma) +
  scale_color_viridis_d(option = "D", begin = 0.3, end = 0.8)
## note that our 7 negative controls all have moderate to high numbers of reads per sample

## create ASV matrix:
asv_mat <- dcast(df, SampleID ~ ASVid, value.var = "Reads", fill=0)
row.names(asv_mat) <- asv_mat$SampleID
asv_mat$SampleID <- NULL

## import metadata and ASV tables into phyloseq
phy_mat <- otu_table(asv_mat, taxa_are_rows = FALSE)
phy_meta <- sample_data(sumry0)
ps <- phyloseq(phy_mat, phy_meta)

## generate stats from decontam
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)  ## ~ 1% of ASVs are flagged as contaminants... not very worrisome

## gather data to make plot of prevalence of TRUE/FALSE contam calls:
contamdf.prev$ASVid <- row.names(contamdf.prev)
ContamSumry <- contamdf.prev %>% select(ASVid, contaminant, prev)

ASVsumry <- df %>% 
  group_by(ASVid, SampleType) %>% 
  summarise(counts=n()) %>%
  spread(ASVsumry, key=SampleType, value=counts, convert = )

contam_df <- merge(ASVsumry, ContamSumry)

ggplot(contam_df, aes(x=prev, y=prev, color=contaminant)) + geom_point()




sumry0$is.neg <- sumry0$SampleType == "control"
contamdf.prev <- isContaminant(asv_mat, method="prevalence", neg="is.neg")


ps.pa.neg <- sumry0 %>% filter(SampleType == "control")

ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
