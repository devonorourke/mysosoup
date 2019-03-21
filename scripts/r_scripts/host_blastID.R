library(taxize)
library(readr)
library(dplyr)

mangan_blast_in <- read_csv(file="https://github.com/devonorourke/mysosoup/raw/master/data/host/mangan_hittable.csv", col_names = FALSE)
#runfromlocal: mangan_blast_in <- read_csv("~/Repos/mysosoup/data/host/mangan_hittable.csv", col_names = FALSE)
colnames(mangan_blast_in) <- c("query", "subject", "pid", "alnlen", "mismatch", "gapopen", 
                               "qstart", "qend", "sstart", "send", "eval", "bitscore")

## get unique IDs to search
mangan_hit_taxIDs <- unique(mangan_blast_in$subject)

## get  NCBI IDs to search
mangan.gb2uid <- genbank2uid(id = mangan_hit_taxIDs, verbose = TRUE)

## collapse list, then manually search the remaining unique records with 'taxize' call...
## passing in vector results in an NCBI timeout error
test <- as.data.frame(do.call("rbind", mangan.gb2uid))
test$V1 <- as.character(test$V1)
searchthese <- unique(test$V1)
searchthese

out1 <- classification(385036, db='ncbi')   ## Myotis sodalis
out2 <- classification(59463, db='ncbi')    ## Myotis lucifugus
out3 <- classification(27670, db='ncbi')    ## Nycticeius humeralis
out4 <- classification(52099, db='ncbi')    ## Desmognathus conanti (salamander!?)    --- very poor alignment identity
out5 <- classification(1656509, db='ncbi')  ## Hyaliodes vitripennis (true bug)
out6 <- classification(159327, db='ncbi')   ## Myotis keaysi (central/south american bat)     --- very poor alignment identity
out7 <- classification(109479, db='ncbi')   ## Myotis mystacinus (european bat)     --- very poor alignment identity

## include results with original data.frame:
test$V2 <- mangan_hit_taxIDs
colnames(test) <- c("ncbiID", "subject")
mangan_blast_out <- merge(mangan_blast_in, test, all.x=TRUE)

mangan_blast_out$taxa <- ""
mangan_blast_out$taxa <- ifelse(grepl('385036', mangan_blast_out$ncbiID), gsub("", "Myotis sodalis", mangan_blast_out$taxa), mangan_blast_out$taxa)
mangan_blast_out$taxa <- ifelse(grepl('59463', mangan_blast_out$ncbiID), gsub("", "Myotis lucifugus", mangan_blast_out$taxa), mangan_blast_out$taxa)
mangan_blast_out$taxa <- ifelse(grepl('27670', mangan_blast_out$ncbiID), gsub("", "Nycticeius humeralis", mangan_blast_out$taxa), mangan_blast_out$taxa)
mangan_blast_out$taxa <- ifelse(grepl('52099', mangan_blast_out$ncbiID), gsub("", "Desmognathus conanti", mangan_blast_out$taxa), mangan_blast_out$taxa)
mangan_blast_out$taxa <- ifelse(grepl('1656509', mangan_blast_out$ncbiID), gsub("", "Hyaliodes vitripennis", mangan_blast_out$taxa), mangan_blast_out$taxa)
mangan_blast_out$taxa <- ifelse(grepl('159327', mangan_blast_out$ncbiID), gsub("", "Myotis keaysi", mangan_blast_out$taxa), mangan_blast_out$taxa)
mangan_blast_out$taxa <- ifelse(grepl('109479', mangan_blast_out$ncbiID), gsub("", "Myotis mystacinus", mangan_blast_out$taxa), mangan_blast_out$taxa)

write.csv(mangan_blast_out, file="~/Repos/mysosoup/data/host/Mangan.blastOut.csv", quote=FALSE, row.names = FALSE)


## create a list of HashIDs from three bat species we believe are hosts:
host_species <- c('Myotis sodalis', 'Myotis lucifugus', 'Nycticeius humeralis')
host_out <- mangan_blast_out %>% filter(taxa %in% host_species) %>% select(query) %>% distinct(.)
colnames(host_out) <- '#OTUID'
write.table(host_out, file = "~/Repos/mysosoup/data/host/host_hashIDs.txt", row.names = FALSE, col.names = TRUE, quote=FALSE)

## quick visualization to confirm that each ASV ID has only one taxa assigned to it
  ## it's possible that these taxa have percent identities similar enough to overlap...
host_tmp <- mangan_blast_out %>% filter(taxa %in% host_species)
## nope - good!