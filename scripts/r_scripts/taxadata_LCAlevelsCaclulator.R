library(tidyverse)

## import data
p97c94 <- read_delim("~/Repos/mysosoup/data/taxonomy/testtax/mangan_arthtax_p97c94.tsv", col_names = TRUE, delim = "\t")
p94c94 <- read_delim("~/Repos/mysosoup/data/taxonomy/testtax/mangan_arthtax_p94c94.tsv", col_names = TRUE, delim = "\t")
p92c90 <- read_delim("~/Repos/mysosoup/data/taxonomy/testtax/mangan_arthtax_p92c90.tsv", col_names = TRUE, delim = "\t")
nb <- read_delim("~/Repos/mysosoup/data/taxonomy/testtax/mangan_arthtax_skl.tsv", col_names = TRUE, delim = "\t")

## how many 'Ambiguous_taxa' are there, per level?
## long function for calculation
LCAfinder.function <- function(data) {
  tmp <- data %>%
    filter(., grepl("Ambiguous_taxa", Taxon)) %>%
    separate(data=., col=Taxon, sep = ";", 
             into = c("kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")) %>%
    gather(Confidence, kingdom_name, phylum_name, class_name, order_name, family_name, genus_name, species_name, value = 'rank') %>%
    mutate(AmbigBool = ifelse(grepl("Ambiguous_taxa", rank), 1, 0))
  
  kingdomAmbig <- with(tmp, tmp[ Confidence == 'kingdom_name' & AmbigBool == 1,])
  phylumAmbig <- with(tmp, tmp[ Confidence == 'phylum_name' & AmbigBool == 1,])
  classAmbig <- with(tmp, tmp[ Confidence == 'class_name' & AmbigBool == 1,])
  orderAmbig <- with(tmp, tmp[ Confidence == 'order_name' & AmbigBool == 1,])
  familyAmbig <- with(tmp, tmp[ Confidence == 'family_name' & AmbigBool == 1,])
  genusAmbig <- with(tmp, tmp[ Confidence == 'genus_name' & AmbigBool == 1,])
  speciesAmbig <- with(tmp, tmp[ Confidence == 'species_name' & AmbigBool == 1,])
  
  speciesLCA <- data.frame(setdiff(speciesAmbig$`Feature ID`, genusAmbig$`Feature ID`)) %>% mutate(LCAlevel="species")
  colnames(speciesLCA)[1] <- "ASVid"
  genusLCA <- data.frame(setdiff(genusAmbig$`Feature ID`, familyAmbig$`Feature ID`)) %>% mutate(LCAlevel="genus")
  colnames(genusLCA)[1] <- "ASVid"
  familyLCA <- data.frame(setdiff(familyAmbig$`Feature ID`, orderAmbig$`Feature ID`)) %>% mutate(LCAlevel="family")
  colnames(familyLCA)[1] <- "ASVid"
  orderLCA <- data.frame(setdiff(orderAmbig$`Feature ID`, classAmbig$`Feature ID`)) %>% mutate(LCAlevel="order")
  colnames(orderLCA)[1] <- "ASVid"
  classLCA <- data.frame(setdiff(classAmbig$`Feature ID`, phylumAmbig$`Feature ID`)) %>% mutate(LCAlevel="class")
  colnames(classLCA)[1] <- "ASVid"
  phylumLCA <- data.frame(setdiff(phylumAmbig$`Feature ID`, kingdomAmbig$`Feature ID`)) %>% mutate(LCAlevel="phylum")
  colnames(phylumLCA)[1] <- "ASVid"
  
  rbind(speciesLCA, genusLCA, familyLCA, orderLCA, classLCA, phylumLCA)
}

p97c94_lcaASVs <- LCAfinder.function(p97c94)  ## just 58 ASVs are Ambiguous (24 Genus, 34 Species)
p94c94_lcaASVs <- LCAfinder.function(p94c94)  ## 6 ASVs with LCA applied (1 Genus, 5 Species)
p92c90_lcaASVs <- LCAfinder.function(p92c90)  ## 3 ASVs with LCA applied (3 Species)
nb_lcaASVs <- LCAfinder.function(nb)          ## 56 ASVs with LCA applied (5 Genus, 49 Species)
