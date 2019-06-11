## !! CAUTION
## ahead of this script, I've added the following line to the .Renviron
## --> ENTREZ_KEY='4868bb1d84231abb79ca9fc41a9cd18ca109'
## That API key was generated from my NCBI account. See https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
## without that key being set in the .Renviron, the process of querying NCBI will fail (it aborts if there are more than 3 request)

library(taxize)
library(tidyverse)

mangan_blast_in <- read_csv(file="https://github.com/devonorourke/mysosoup/raw/master/data/host/mangan_hittable.csv.gz", col_names = FALSE)
#runfromlocal: mangan_blast_in <- read_csv("~/Repos/mysosoup/data/host/mangan_hittable2.csv.gz", col_names = FALSE)
colnames(mangan_blast_in) <- c("query", "subject", "pid", "alnlen", "mismatch", "gapopen", 
                               "qstart", "qend", "sstart", "send", "eval", "bitscore")

## get unique IDs to search
mangan_hit_taxIDs <- unique(mangan_blast_in$subject)

## get  NCBI IDs to search
mangan.gb2uid <- genbank2uid(id = mangan_hit_taxIDs, verbose = TRUE)

## gather NCBI taxa
test <- as.data.frame(do.call("rbind", mangan.gb2uid))  ## the NCBI accession values
test$V1 <- as.character(test$V1)
searchthese <- unique(test$V1)                          ## the NCBI tax IDs
classify_list <- classification(searchthese, db='ncbi') ## the NCBI taxa names

classify_df <- do.call(rbind.data.frame, classify_list) ## list to data.frame
classify_df$tmp <- row.names(classify_df)
rankFilter <- c("family", "genus", "species")
classify_out <- classify_df %>% 
  separate(data = ., col = 'tmp', into = c('ncbiID', 'drop'), sep = '\\.') %>% 
  select(-drop, -id) %>% 
  filter(rank %in% rankFilter) %>% 
  spread(rank, name)                                    ## spread into family-genus-species info, per Subject

## include results with original data.frame:
test$V2 <- mangan_hit_taxIDs
colnames(test) <- c("ncbiID", "subject")
mangan_blast_out <- merge(mangan_blast_in, test, all.x=TRUE)
mangan_blast_out <- merge(mangan_blast_out, classify_out)

write.csv(mangan_blast_out, file="~/Repos/mysosoup/data/host/Mangan.blastOut.csv", quote=FALSE, row.names = FALSE)


## create a list of HashIDs from three bat species we believe are hosts:
host_species <- c('Myotis sodalis', 'Myotis lucifugus', 'Nycticeius humeralis')
host_out <- mangan_blast_out %>% filter(species %in% host_species) %>% select(query) %>% distinct(.)
colnames(host_out) <- '#OTUID'
write.table(host_out, file = "~/Repos/mysosoup/data/host/host_hashIDs_ncbi.txt", row.names = FALSE, col.names = TRUE, quote=FALSE)

## quick visualization to confirm that each ASV ID has only one taxa assigned to it
  ## it's possible that these taxa have percent identities similar enough to overlap...
host_tmp <- mangan_blast_out %>% filter(species %in% host_species)
ggplot(host_tmp, aes(x=query, y=species)) + geom_point() + theme(axis.text.x = element_blank())
  ## nope - good!
