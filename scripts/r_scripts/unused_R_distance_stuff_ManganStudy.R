## Try ordinating with PCoA instead of NMDS using APE package:
pcoaFracVarfunction <- function(distance,DistMethod) {
  pcoa_list <- pcoa(distance)
  pcoa_Eigvals <- pcoa_list$values
  pcoa_EigData <- pcoa_Eigvals %>% select(Cum_corr_eig) %>%
    mutate(Axis=paste0("Axis.",row.names(.))) %>%
    mutate(tmpVal=c(.$Cum_corr_eig[1],diff(.$Cum_corr_eig))) %>%
    mutate(FracVar=round(tmpVal*100, 2)) %>%
    select(FracVar, Axis) %>%
    mutate(DistMethod=DistMethod)
}

pcoaPointsfunction <- function(distance,DistMethod) {
  pcoa_list <- pcoa(distance)
  pcoa_Points <- data.frame(pcoa_list$vectors) %>%
    mutate(SampleID = row.names(.)) %>%
    mutate(DistMethod = DistMethod)
}


raup.pcoa.fracvar <- pcoaFracVarfunction(raupDist, "Raup")
raup.pcoa.points <- pcoaPointsfunction(raupDist, "Raup")
raup.pcoa.plot <- merge(raup.pcoa.points, tinymeta)
AxisX <- raup.pcoa.fracvar %>% filter(Axis == "Axis.1") %>% select(FracVar) %>% pull()
AxisY <- raup.pcoa.fracvar %>% filter(Axis == "Axis.2") %>% select(FracVar) %>% pull()

ggplot() +
  geom_point(data = raup.pcoa.plot, aes(x=Axis.1, y=Axis.2, color=CollectionMonth)) +
  labs(x=paste("Axis.1 [", AxisX, "]"),
       y=paste("Axis.2 [", AxisY, "]"))



bray.pcoa.fracvar <- pcoaFracVarfunction(brayDist, "Bray")
bray.pcoa.points <- pcoaPointsfunction(brayDist, "Bray")
bray.pcoa.plot <- merge(bray.pcoa.points, tinymeta)
AxisX <- bray.pcoa.fracvar %>% filter(Axis == "Axis.1") %>% select(FracVar) %>% pull()
AxisY <- bray.pcoa.fracvar %>% filter(Axis == "Axis.2") %>% select(FracVar) %>% pull()

ggplot() +
  geom_point(data = bray.pcoa.plot, aes(x=Axis.1, y=Axis.2, color=CollectionMonth)) +
  labs(x=paste("Axis.1 [", AxisX, "]"),
       y=paste("Axis.2 [", AxisY, "]"))



mori.pcoa.fracvar <- pcoaFracVarfunction(moriDist, "Morisita")
mori.pcoa.points <- pcoaPointsfunction(moriDist, "Morisita")
mori.pcoa.plot <- merge(mori.pcoa.points, tinymeta)
AxisX <- mori.pcoa.fracvar %>% filter(Axis == "Axis.1") %>% select(FracVar) %>% pull()
AxisY <- mori.pcoa.fracvar %>% filter(Axis == "Axis.2") %>% select(FracVar) %>% pull()

ggplot() +
  geom_point(data = mori.pcoa.plot, aes(x=Axis.1, y=Axis.2, color=CollectionMonth)) +
  labs(x=paste("Axis.1 [", AxisX, "]"),
       y=paste("Axis.2 [", AxisY, "]"))




plot_scree(ordinate(mr, "PCoA", "bray"))
plot_scree(ordinate(mr, "NMDS", "bray"))
plot_scree(ordinate(mr, "DCA", "bray", ))

?ordinate()


## Estimate distances in Phyloseq with DCA method
p.raupDCA  <- ordinate(mr, distance=raupDist)
p.brayDCA  <- ordinate(mr, distance=brayDist)
p.moriDCA  <- ordinate(mr, method = "DCA", distance=moriDist)

p.raup_df <- data.frame(p.raupDCA$rproj) %>% mutate(SampleID=row.names(.)) %>% mutate(Dist="raup")
p.bray_df <- data.frame(p.brayDCA$rproj) %>% mutate(SampleID=row.names(.)) %>% mutate(Dist="bray")
p.mori_df <- data.frame(p.moriDCA$rproj) %>% mutate(SampleID=row.names(.)) %>% mutate(Dist="mori")
p.DCA.df <- rbind(p.raup_df, p.bray_df, p.mori_df)
rm(p.raup_df, p.bray_df, p.mori_df)
p.DCA.df <- merge(p.DCA.df, meta)

## plot
ggplot(p.DCA.df, aes(x=DCA1, y=DCA2, label=SampleID)) + geom_point() + facet_grid(~Dist)

farRight <- '9152017HBD4'
topLeft <- '9152017EGB9'
botLeft <- '7272017HBA2'
DCAoutlierSamples <- c(farRight, topLeft, botLeft)

data_outliers <- famData %>% filter(SampleID %in% DCAoutlierSamples)


library(qiime2R)
library(ape)

uuni_fam <- read_qza("~/Repos/mysosoup/data/qiime_qza/distances/mangan_nobats_uuni_dist_FamOnly.qza")
uuni_fam_dist <- uuni_fam$data
wuni_fam <- read_qza("~/Repos/mysosoup/data/qiime_qza/distances/mangan_nobats_wuni_dist_FamOnly.qza")
wuni_fam_dist <- wuni_fam$data
rm(uuni_fam, wuni_fam)

pcoaFracVarfunction <- function(distance,DistMethod) {
  pcoa_list <- pcoa(distance)
  pcoa_Eigvals <- pcoa_list$values
  pcoa_EigData <- pcoa_Eigvals %>% select(Cum_corr_eig) %>%
    mutate(Axis=paste0("Axis.",row.names(.))) %>%
    mutate(tmpVal=c(.$Cum_corr_eig[1],diff(.$Cum_corr_eig))) %>%
    mutate(FracVar=round(tmpVal*100, 2)) %>%
    select(FracVar, Axis) %>%
    mutate(DistMethod=DistMethod)
}

pcoaPointsfunction <- function(distance,DistMethod) {
  pcoa_list <- pcoa(distance)
  pcoa_Points <- data.frame(pcoa_list$vectors) %>%
    mutate(SampleID = row.names(.)) %>%
    mutate(DistMethod = DistMethod)
}


uuni.pcoa.fracvar <- pcoaFracVarfunction(uuni_fam_dist, "Uuni")
uuni.pcoa.points <- pcoaPointsfunction(uuni_fam_dist, "Uuni")
uuni.pcoa.points <- merge(uuni.pcoa.points, tinymeta)
AxisX <- uuni.pcoa.fracvar %>% filter(Axis == "Axis.1") %>% select(FracVar) %>% pull()
AxisY <- uuni.pcoa.fracvar %>% filter(Axis == "Axis.2") %>% select(FracVar) %>% pull()

ggplot() +
  geom_point(data = uuni.pcoa.points, aes(x=Axis.1, y=Axis.2, color=Site, shape=CollectionMonth)) +
  labs(x=paste("Axis.1 [", AxisX, "]"),
       y=paste("Axis.2 [", AxisY, "]"))


wuni.pcoa.fracvar <- pcoaFracVarfunction(wuni_fam_dist, "Wuni")
wuni.pcoa.points <- pcoaPointsfunction(wuni_fam_dist, "Wuni")
wuni.pcoa.points <- merge(wuni.pcoa.points, tinymeta)
AxisX <- wuni.pcoa.fracvar %>% filter(Axis == "Axis.1") %>% select(FracVar) %>% pull()
AxisY <- wuni.pcoa.fracvar %>% filter(Axis == "Axis.2") %>% select(FracVar) %>% pull()

ggplot() +
  geom_point(data = wuni.pcoa.points, aes(x=Axis.1, y=Axis.2, color=CollectionMonth)) +
  labs(x=paste("Axis.1 [", AxisX, "]"),
       y=paste("Axis.2 [", AxisY, "]"))
