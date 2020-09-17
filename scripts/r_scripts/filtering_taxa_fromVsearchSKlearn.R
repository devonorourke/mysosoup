library(tidyverse)

################################################################################
## step 1 - Import all VSEARCH classified clustered seqs and filter candidate features
################################################################################
## starting with 1,936 sequence clusters (98.5% identity from denoised ASVs)

## import vsearch alignments that required 100% identity across >= 94% query coverage
## LCA process applied if sequence cluster met criteria threshold
vs_raw <- read_delim(file="https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/taxonomy/vsearchOnly_p100_c94_taxonomy.tsv",
                     delim = "\t", col_names = TRUE)

## remove 'Ambiguous taxa' from strings
vs_raw$Taxon <- gsub("Ambiguous_taxa", "", vs_raw$Taxon)

## remove sequences without any assignment
vs_filt <- vs_raw %>% 
  filter(Taxon != "Unassigned")
  ## 556 remaining seqs

## split 'Taxon' field into kingdom-->species levels
vs_filt <- vs_filt %>% 
  separate(Taxon,
           into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

## retain only sequences in phylum Arthropoda
vs_filt <- vs_filt %>% 
  filter(Phylum=="Arthropoda")
  ## 550 seqs remain

## retain only sequences with at least Family-rank information
vs_filt <- vs_filt %>% 
  mutate_all(na_if,"") %>% 
  filter(!is.na(Order)) %>% 
  filter(!is.na(Family))
  ## 549 seqs remain
  ## These 549 seqs represent our VSEARCH candidates

## add classifier label
vs_filt$Classifier <- "VSEARCH"

## create list of FeatureIDs to filter from qiime .qza taxonomy artifact
vslist <- as.data.frame(vs_filt$`Feature ID`)
colnames(vslist) <- 'Feature ID'
write.csv(vslist, file="~/github/mysosoup/data/taxonomy/filtd_taxlist_vsearch.txt",
          row.names = FALSE, quote = FALSE)

## cleanup
rm(vs_raw, vslist)

################################################################################
## step 2 - Import all sklearn classified clustered, ...
## ...then gather only those seqs not retained by VSEARCH above ...
## ...and filter using the same principles applied above
################################################################################

## import sklearn classified seq features
sk_raw <- read_delim(file="~/github/mysosoup/data/taxonomy/sklearn_taxonomy.tsv",
                     delim = "\t", col_names = TRUE)

## exclude vsearch features retained above from sklearn
sk_filt <- sk_raw %>% 
  filter(!`Feature ID` %in% vs_filt$`Feature ID`)
  ## 1387 features remain

## filter using same approach as in VSEARCH above
sk_filt$Taxon <- gsub("Ambiguous_taxa", "", sk_filt$Taxon)

sk_filt <- sk_filt %>% 
  filter(Taxon != "Unassigned")
## 1264 features remain

sk_filt <- sk_filt %>% 
  separate(Taxon,
           into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

sk_filt <- sk_filt %>% 
  filter(Phylum=="Arthropoda")
  ## 1021 seqs remain

sk_filt <- sk_filt %>% 
  mutate_all(na_if,"") %>% 
  filter(!is.na(Order)) %>% 
  filter(!is.na(Family))
  ## 575 features remain

## add classifier label
sk_filt$Classifier <- "sklearn"

## create list of FeatureIDs to filter from qiime .qza taxonomy artifact
sklist <- as.data.frame(sk_filt$`Feature ID`)
colnames(sklist) <- 'Feature ID'
write.csv(sklist, file="~/github/mysosoup/data/taxonomy/filtd_taxlist_sklearn.txt",
          row.names = FALSE, quote = FALSE)

## cleanup
rm(sk_raw, sklist)

################################################################################
## step 3 - combine taxonomy results to generate data table
## also, perform sanity check to ensure featureIDs are distinct in both sklearn and vsearch
################################################################################

## combine vsearch and sklearn into single data frame
all_dat <- bind_rows(sk_filt, vs_filt)

## sanity check we don't have duplicate Feature IDs
finalrows <- nrow(all_dat)
uniquerows <- length(unique(all_dat$`Feature ID`))
finalrows == uniquerows  ## all good!

## write to disk
write.csv(all_dat, file="~/github/mysosoup/data/taxonomy/filtd_tax_dataframe_ALL.csv",
          row.names = FALSE, quote = FALSE)
