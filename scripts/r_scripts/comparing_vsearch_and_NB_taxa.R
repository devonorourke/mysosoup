library(tidyverse)

## import data
p97c94 <- read_delim("~/Repos/mysosoup/data/taxonomy/testtax/mangan_arthtax_p97c94.tsv", col_names = TRUE, delim = "\t")
p94c94 <- read_delim("~/Repos/mysosoup/data/taxonomy/testtax/mangan_arthtax_p94c94.tsv", col_names = TRUE, delim = "\t")
p92c90 <- read_delim("~/Repos/mysosoup/data/taxonomy/testtax/mangan_arthtax_p92c90.tsv", col_names = TRUE, delim = "\t")
nb <- read_delim("~/Repos/mysosoup/data/taxonomy/testtax/mangan_arthtax_skl.tsv", col_names = TRUE, delim = "\t")

## cleanup function
cleanup.function <- function(data) {
  tmp <- data %>% 
    separate(., col = Taxon, sep=';', 
             into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"))
  tmp <- as.data.frame(apply(tmp, 2, function(y) gsub(".__", "", y)))
  tmp <- as.data.frame(apply(tmp, 2, function(y) gsub("^$|^ $", NA, y)))
  tmp <- as.data.frame(apply(tmp, 2, function(y) gsub("Ambiguous_taxa", NA, y)))
}

clean_noAmbig_p97c94 <- cleanup.function(p97c94)
clean_noAmbig_p94c94 <- cleanup.function(p94c94)
clean_noAmbig_p92c90 <- cleanup.function(p92c90)
clean_noAmbig_nb <- cleanup.function(nb)

## retain datasets with at least Family-rank information:
famonly_p97c94 <- clean_noAmbig_p97c94 %>% filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
  filter(!is.na(order_name)) %>% filter(!is.na(family_name))  ## preserves 2990 asvs (about 56.5%)
famonly_nb <- clean_noAmbig_nb %>% filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) %>%
  filter(!is.na(order_name)) %>% filter(!is.na(family_name))  ## preserves 3049 asvs (about 57.6%)

## what's different/common between these datasets?
length(intersect(famonly_nb$`Feature ID`, famonly_p97c94$`Feature ID`)) ## share 2477 ASVs in common (so about 16% are missing from a given set)
length(setdiff(famonly_nb$`Feature ID`, famonly_p97c94$`Feature ID`)) ## 572 missing from p97c94
length(setdiff(famonly_p97c94$`Feature ID`, famonly_nb$`Feature ID`)) ## 513 missing from nb


## retain datasets with at least Class-rank information:
classonly_p97c94 <- clean_noAmbig_p97c94 %>% filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) ## 3082 ASVs (retains 58.3%)
classonly_nb <- clean_noAmbig_nb %>% filter(!is.na(phylum_name)) %>% filter(!is.na(class_name)) ## 3082 ASVs (retains 92.9%)
  ## big distinction: keeping at Class preserves almost all ASVs for Naive Bayes classifier

## what's diff/comm here?
length(intersect(classonly_nb$`Feature ID`, classonly_p97c94$`Feature ID`)) ## share 3080 ASVs in common
length(setdiff(classonly_nb$`Feature ID`, classonly_p97c94$`Feature ID`)) ## 1835 missing from p97c94
length(setdiff(classonly_p97c94$`Feature ID`, classonly_nb$`Feature ID`)) ## just 2 missing from nb

## --------------------------------------------  ##

## What are the read depths and frequencies of ASVs in these Class/Family-filtered datasets (for those in and out of the dataset)?
data <- read_csv(file = "https://github.com/devonorourke/mysosoup/raw/master/data/mangan.asvtable.long.csv.gz", col_names = TRUE)
asvSumry <- data %>% group_by(ASVid) %>% summarise(ReadSums=sum(Reads), nSamples=n())

ggplot(asvSumry, aes(x=ReadSums, y=nSamples)) +
  geom_point(data = asvSumry %>% filter(ASVid %in% classonly_p97c94$`Feature ID`), color="red") +
  geom_point(data = asvSumry %>% filter(!ASVid %in% classonly_p97c94$`Feature ID`), color="blue")

ggplot(asvSumry, aes(x=ReadSums, y=nSamples)) +
  geom_point(data = asvSumry %>% filter(ASVid %in% famonly_p97c94$`Feature ID`), color="red") +
  geom_point(data = asvSumry %>% filter(!ASVid %in% famonly_p97c94$`Feature ID`), color="blue")

ggplot(asvSumry, aes(x=ReadSums, y=nSamples)) +
  geom_point(data = asvSumry %>% filter(ASVid %in% classonly_nb$`Feature ID`), color="black", alpha=0.5) +
  geom_point(data = asvSumry %>% filter(!ASVid %in% classonly_nb$`Feature ID`), color="blue")

ggplot(asvSumry, aes(x=ReadSums, y=nSamples)) +
  geom_point(data = asvSumry %>% filter(ASVid %in% famonly_nb$`Feature ID`), color="red") +
  geom_point(data = asvSumry %>% filter(!ASVid %in% famonly_nb$`Feature ID`), color="blue")

### --------------------------

classonly_p97c94_ASVsumry <- data %>% 
  group_by(ASVid) %>% 
  filter(!ASVid %in% classonly_p97c94$`Feature ID`) %>%
  summarise(meanASVread=mean(Reads), medianASVread=median(Reads), sumASVreads=sum(Reads), nSamplesASV=n())

classonly_nb_ASVsumry <- data %>% 
  group_by(ASVid) %>% 
  filter(!ASVid %in% classonly_nb$`Feature ID`) %>%
  summarise(meanASVread=mean(Reads), medianASVread=median(Reads), sumASVreads=sum(Reads), nSamplesASV=n())
