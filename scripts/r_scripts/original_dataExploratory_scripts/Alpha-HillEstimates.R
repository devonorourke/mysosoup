library(tidyverse)
library(vegan)
library(scales)
library(qiime2R)
library(reshape2)
library(formattable)
library(ggpubr)
library(Matrix)
library(multcompView)

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
meta <- meta %>% select(SampleID, Roost, CollectionMonth, Site, SampleType, BatchType)
meta$Site <- ifelse(meta$Site == "Egner", gsub("Egner", "EN", meta$Site), meta$Site)
meta$Site <- ifelse(meta$Site == "HickoryBottoms", gsub("HickoryBottoms", "HB", meta$Site), meta$Site)
meta$CollectionMonth[is.na(meta$CollectionMonth)] <- "control"
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "6", gsub("6", "June", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "7", gsub("7", "July", meta$CollectionMonth), meta$CollectionMonth)
meta$CollectionMonth <- ifelse(meta$CollectionMonth == "9", gsub("9", "September", meta$CollectionMonth), meta$CollectionMonth)
meta$Labeler <- paste(meta$Site, meta$Roost, sep="-")
meta$Labeler <- ifelse(meta$Labeler == "control-control", gsub("control-control", "control", meta$Labeler), meta$Labeler)
meta$CollectionMonth <- as.factor(meta$CollectionMonth)


## add taxonomy information
taxa <- read_delim(file = "https://github.com/devonorourke/mysosoup/raw/master/data/taxonomy/mangan_tax_p97c94.tsv", delim = "\t", col_names = TRUE)
taxa <- taxa %>% separate(., col = Taxon, sep=';', into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>% select(-Confidence)
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub(".__", "", y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("^$|^ $", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
taxa <- as.data.frame(apply(taxa, 2, function(y) gsub("Unassigned", NA, y)))
colnames(taxa)[1] <- "ASVid"

## import dataset, filter out reads all non-Arthropod reads without (at least) Family-rank info, ..
# alternatively: download files:
## download.file("https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza", "tmp.qza")
# then set PATH to whatever directory needed:
# qzapath = "PATH/TO/tmp.qza"

## ..filter out remaining ASVs present exclusively in control samples, then calculate Hill Number values
## run from local: qzapath = "~/Repos/mysosoup/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza"
featuretable <- read_qza(qzapath)
mat.tmp <- featuretable$data
rm(featuretable)
df.tmp <- as.data.frame(mat.tmp)
rm(mat.tmp)
df.tmp$OTUid <- rownames(df.tmp)
rownames(df.tmp) <- NULL
tmp <- melt(df.tmp, id = "OTUid") %>% filter(value != 0)
rm(df.tmp) 
colnames(tmp) <- c("ASVid", "SampleID", "Reads")
df.tmp <- merge(tmp, taxa)
rm(tmp)
df.tmp <- merge(df.tmp, meta)
onlyControlASVs <- setdiff(df.tmp %>% filter(SampleType=="control") %>% select(ASVid) %>% pull(), df.tmp %>% filter(SampleType=="sample") %>% select(ASVid) %>% pull())
df_filt.tmp <- df.tmp %>% filter(!ASVid %in% onlyControlASVs) %>% filter(phylum_name=="Arthropoda") %>% filter(!is.na(family_name)) %>% filter(SampleType == "sample")
rm(df.tmp)
mat.tmp <- dcast(data = df_filt.tmp, formula = SampleID ~ ASVid, value.var='Reads', fill = 0)
row.names(mat.tmp) <- mat.tmp$SampleID
mat.tmp$SampleID <- NULL
tmp.hill <- data.frame(renyi(mat.tmp, scales = c(0,1,2), hill=TRUE)) %>% mutate(SampleID = row.names(.))
colnames(tmp.hill)[1:3] <- c("q=0", "q=1", "q=2")
Alpha_df <- gather(tmp.hill, key="Hill_qType", value = "Hill_value", c('q=0', 'q=1', 'q=2'))
mangan_hill_df <- merge(Alpha_df, meta)
rm(tmp.hill, mat.tmp, Alpha_df)


## plot parameter setup
mangan_hill_df$CollectionMonth <- factor(mangan_hill_df$CollectionMonth, levels = c("June", "July", "September"))
v3pal <- viridis::plasma(3, begin = 0.35, end = 0.9, direction = -1)

mangan_hill_df$Label2 <- paste(mangan_hill_df$Site, mangan_hill_df$CollectionMonth, sep="-")

## plot alpha diversity collectively
## save as all_Alpha_Hillvals; export at 900x697
ggplot(mangan_hill_df, aes(x=Site, y=Hill_value, color=CollectionMonth)) + 
  geom_boxplot(outlier.shape = NA, color="gray30") +
  geom_jitter(width = 0.2, alpha=0.6) + 
  scale_color_manual(values = c(v3pal, "gray40"), labels=c("June", "July", "September")) +
  facet_grid(Hill_qType ~ CollectionMonth) +
  labs(x="", y="Estimated diversity", color = "Month") +
  theme_devon() + theme(legend.position = "none")

## alternative plot groups all data by Site, not by Roost
## save as all_Alpha_Hillvals_separateRoosts; export at 900x697
ggplot(mangan_hill_df, aes(x=Labeler, y=Hill_value, color=CollectionMonth)) + 
  geom_jitter(width = 0.2, alpha=0.8) + 
  scale_color_manual(values = c(v3pal, "gray40"), labels=c("June", "July", "September")) +
  facet_wrap(~ Hill_qType ) +
  labs(x="", y="Estimated diversity", color = "Month") +
  theme_devon() +
  theme(axis.text.x = element_text(angle = 22.5, hjust=1),
        legend.position = "top")

##########
## collect some mean/sd stats per Site/Month for each Hill Val
##########
hillstats <- mangan_hill_df %>% 
  group_by(Site, CollectionMonth, Hill_qType) %>% 
  summarise(meanQ=mean(Hill_value), sdQ=sd(Hill_value)) %>% 
  mutate(meanQ=round(meanQ, 2), sdQ=round(sdQ, 2)) %>% 
  rename(Month=CollectionMonth) %>% 
  arrange(Hill_qType, Site, 
          factor(Month, levels = c("June", "July", "September")))
  #arrange(Hill_qType, Site, Month)
  
write.csv(hillstats, file="~/Repos/mysosoup/data/text_tables/Hill_mean_and_SD_stats.csv",
          row.names = FALSE, quote = FALSE)
## or make pretty table; save as 'Hill_mean_and_SD_stats' in figure; export at 450x625
formattable(hillstats)


##### Run ANOVA for group significance:
#### Run ANOVA for group significance:
anovafunction <- function(data, qval, filename) {
  aov.out <- aov(Hill_value ~ CollectionMonth * Site, data=data %>% filter(Hill_qType == qval))
  capture.output(summary(aov.out),file=paste0("~/Repos/mysosoup/data/text_tables/anovas/",filename))
}

anovafunction(mangan_hill_df, "q=0", "aov_mangan_alpha_q0.txt")
anovafunction(mangan_hill_df, "q=1", "aov_mangan_alpha_q1.txt")
anovafunction(mangan_hill_df, "q=2", "aov_mangan_alpha_q2.txt")


###### Run Kruskal Wallis for nonparametric test of same data:
kwfunction <- function(data, qval){
  data$Grouper <- paste(data$CollectionMonth, data$Site, sep="-")
  data$Grouper <- as.factor(data$Grouper)
  kruskal.test(Hill_value ~ Grouper, data = data %>% filter(Hill_qType == qval))
}

capture.output(kwfunction(mangan_hill_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q0.txt")
capture.output(kwfunction(mangan_hill_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q1.txt")
capture.output(kwfunction(mangan_hill_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/kruskal/kw_mangan_alpha_q2.txt")


## Dunn test 
library(FSA)

dunnfunction <- function(data, qfilt){
  data$Grouper <- paste(data$CollectionMonth, data$Site, sep="-")
  data$Grouper <- as.factor(data$Grouper)
  tmp <- dunnTest(Hill_value ~ Grouper, data=data %>% filter(Hill_qType==qfilt), method = "bh") 
  tmp <- data.frame(tmp$res)
  tmp <- tmp %>% mutate(P.adj=round(P.adj, 3))
  tmp <- tmp %>% mutate(Z=round(Z, 3))
  tmp <- tmp %>% mutate(P.unadj=round(P.unadj, 3))
  tmp %>% arrange(P.adj)
}

capture.output(dunnfunction(mangan_hill_df, "q=0"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q0.txt")
capture.output(dunnfunction(mangan_hill_df, "q=1"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q1.txt")
capture.output(dunnfunction(mangan_hill_df, "q=2"),file="~/Repos/mysosoup/data/text_tables/dunn/dunn_mangan_alpha_q2.txt")


##########
##alternative plot will put significance values over different pairwise comparisons

## split into 3 separate plots; one for each Hill Number
hill0 <- mangan_hill_df %>% filter(Hill_qType=="q=0")
hill0$Label2 <- factor(hill0$Label2, levels = c(
  "EN-June", "HB-June", "EN-July", "HB-July", "EN-September", "HB-September"))

hill1 <- mangan_hill_df %>% filter(Hill_qType=="q=1")
hill1$Label2 <- factor(hill1$Label2, levels = c(
  "EN-June", "HB-June", "EN-July", "HB-July", "EN-September", "HB-September"))

hill2 <- mangan_hill_df %>% filter(Hill_qType=="q=2")
hill2$Label2 <- factor(hill2$Label2, levels = c(
  "EN-June", "HB-June", "EN-July", "HB-July", "EN-September", "HB-September"))

## ggpubr to create individual figures then stitch together after
p0 <- ggboxplot(hill0, x = "Label2", y = "Hill_value",
                color = "CollectionMonth", palette = v3pal,
                add = "jitter") 

## which vals are significant?
compare_means(Hill_value ~ Label2,  data = hill0) %>% arrange(p.format)
q0_comparisons <- list( c("EN-June","EN-September"), 
                        c("EN-July","HB-September"), 
                        c("HB-June","HB-September"),
                        c("EN-July","EN-September"), 
                        c("HB-July","HB-September"), 
                        c("HB-June", "EN-September"),
                        c("EN-June","HB-September"))
## full plot
a <- p0 + 
  labs(y= "Estimated diversity") +
  theme_devon() +
  theme(legend.position="none", 
        axis.text.x = element_text(angle=22.5, hjust=1),
        axis.title.x = element_blank()) +
  #stat_compare_means(comparisons = q0_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 250)     # Add global p-value

################################################################################
## repeat for q=1

p1 <- ggboxplot(hill1, x = "Label2", y = "Hill_value",
                color = "CollectionMonth", palette = v3pal,
                add = "jitter") 

compare_means(Hill_value ~ Label2,  data = hill1) %>% arrange(p.format)

q1_comparisons <- list( c("HB-July","EN-September"), 
                        c("HB-July","HB-September"), 
                        c("EN-June","HB-July")) 

b <- p1 + 
  labs(y= "Estimated diversity") +
  theme_devon() +
  theme(legend.position="none", 
        axis.text.x = element_text(angle=22.5, hjust=1),
        axis.title.x = element_blank()) +
  stat_compare_means(comparisons = q1_comparisons) +
  stat_compare_means(label.y = 70)


################################################################################
## repeat for q=2

p2 <- ggboxplot(hill2, x = "Label2", y = "Hill_value",
                color = "CollectionMonth", palette = v3pal,
                add = "jitter") 

compare_means(Hill_value ~ Label2,  data = hill2) %>% arrange(p.format)

q2_comparisons <- list( c("HB-July","EN-September"), 
                        c("HB-July","HB-September"), 
                        c("EN-June","HB-July"))

c <- p2 + 
  labs(y= "Estimated diversity") +
  theme_devon() +
  theme(legend.position="none", 
        axis.text.x = element_text(angle=22.5, hjust=1),
        axis.title.x = element_blank()) +
  stat_compare_means(comparisons = q2_comparisons) +
  stat_compare_means(label.y = 40)


################################################################################
## stitch all the plots together; save as 'all_Alpha_Hillvals_wWilcox'; export at 550x1100
ggarrange(a,b,c, 
          nrow = 3, ncol=1, 
          #label.x = 0,
          #label.y = 0,
          labels=c("A", "B", "C"))


##########
## Nick wanted to use letters to signify different means instead of the bars and p-values, so we can do that this way:

lpf <- function(data){
  attach(data)
  tmp_tri <- pairwise.wilcox.test(Hill_value, Label2, p.adjust.method = "bonf")
  detach()
  
  h0_mat <- tmp_tri$p.value
  `EN-June`= c(1, NA, NA, NA, NA, NA)
  `HB-September` = c(NA, NA, NA, NA, 1)
  h0_mat <- cbind(h0_mat, `HB-September`)
  h0_mat <- rbind(`EN-June`, h0_mat)
  h0_mat <- Matrix::forceSymmetric(h0_mat,uplo="L")
  h0_mat <- as.matrix(h0_mat)
  h0_mat[is.na(h0_mat)] <- 1
  
  tmp_lmat <- multcompLetters(h0_mat, compare="<=", threshold=0.05, Letters=letters)
  tmp_labels <- data.frame(Letters=tmp_lmat$monospacedLetters) %>% 
    mutate(Label2=row.names(.)) %>% 
    mutate(Letters=trimws(Letters))
}

## Apply function to each Hill Number; may need to adjust the Hill_value (see comment)
hill0_labels <- lpf(hill0)
hill0_labels$Hill_value = 150  ## derive this value from wherever plotting space is available

d <- p0 + 
  stat_compare_means(label.y = 170) +
  geom_text(data = hill0_labels, aes(x=Label2, y=Hill_value, label=Letters)) +
  labs(y= "Estimated diversity") +
  theme_devon() +
  theme(legend.position="none", 
        axis.text.x = element_text(angle=22.5, hjust=1),
        axis.title.x = element_blank())


hill1_labels <- lpf(hill1)
hill1_labels$Hill_value = 53

e <- p1 + 
  stat_compare_means(label.y = 60) +
  geom_text(data = hill1_labels, aes(x=Label2, y=Hill_value, label=Letters)) +
  labs(y= "Estimated diversity") +
  theme_devon() +
  theme(legend.position="none", 
        axis.text.x = element_text(angle=22.5, hjust=1),
        axis.title.x = element_blank())


hill2_labels <- lpf(hill2)
hill2_labels$Hill_value = 30

f <- p2 + 
  stat_compare_means(label.y = 35) +
  geom_text(data = hill2_labels, aes(x=Label2, y=Hill_value, label=Letters)) +
  labs(y= "Estimated diversity") +
  theme_devon() +
  theme(legend.position="none", 
        axis.text.x = element_text(angle=22.5, hjust=1),
        axis.title.x = element_blank())

## stitch all the plots together; save as 'all_Alpha_Hillvals_wLetters'; export at 550x1100
ggarrange(d,e,f, 
          nrow = 3, ncol=1, 
          #label.x = 0,
          #label.y = 0,
          labels=c("A", "B", "C"))
