library(ggmap)
library(ggsn)
library(tidyverse)
library(ggrepel)
library(scico)
library(svglite)

## making map of number of samples sequenced in dataset (regardless of read depth)

sitedata <- read_csv("https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/metadata/brandenbark_pole_coordinates.csv") %>% 
  rename(lat = lat_dec, lon = long_dec) %>% select(-x, -y)

## import map layers as basic stamen map
# sitemap <- get_stamenmap(bbox = c(left = -89.2, bottom = 37.2, right = -88.9, top = 37.4),
#                          zoom = 11,
#                          maptype = "terrain")

############ 

googmap_hyb <- get_googlemap(center=c(lon = -89.075, lat = 37.325), zoom=13, maptype = "hybrid")
googmap_sat <- get_googlemap(center=c(lon = -89.075, lat = 37.325), zoom=13, maptype = "satellite")

ggmap(googmap_sat) +
  geom_point(aes(x=lon, y=lat, color=Area), data=sitedata, size=2.5) +
  scale_color_manual(values=c("gold2", "pink")) +
  annotate("text", x = -89.075, y = 37.345, label = "Hickory", color="pink", size=6) +
  annotate("text", x = -89.1, y = 37.298, label = "Egner", color="gold2", size=6) +
  # geom_label_repel(aes(x=lon, y=lat, label=Pole), data=sitedata %>% filter(Pole=="HB-A"), size=3.5, direction = "x", nudge_x = -0.01) +
  # geom_label_repel(aes(x=lon, y=lat, label=Pole), data=sitedata %>% filter(Pole=="HB-B"), size=3.5, direction = "x", nudge_y = 0.01) +
  # geom_label_repel(aes(x=lon, y=lat, label=Pole), data=sitedata %>% filter(Pole=="HB-C"), size=3.5, direction = "x", nudge_x = 0.01) +
  # geom_label_repel(aes(x=lon, y=lat, label=Pole), data=sitedata %>% filter(Pole=="HB-D"), size=3.5, direction = "x", nudge_y = -0.01) +
  # geom_label_repel(aes(x=lon, y=lat, label=Pole), data=sitedata %>% filter(Pole=="EG-A"), size=3.5, direction = "x", nudge_x = 0.01) +
  # geom_label_repel(aes(x=lon, y=lat, label=Pole), data=sitedata %>% filter(Pole=="EG-B"), size=3.5, direction = "x", nudge_y = 0.01) +
  # geom_label_repel(aes(x=lon, y=lat, label=Pole), data=sitedata %>% filter(Pole=="EG-C"), size=3.5, direction = "x", nudge_y = 0.01) +
  # geom_label_repel(aes(x=lon, y=lat, label=Pole), data=sitedata %>% filter(Pole=="EG-D"), size=3.5, direction = "x", nudge_x = -0.01) +
  scalebar(x.min = -89.05, x.max = -89.075, y.min = 37.285, y.max = 37.29,
           dist = 1, dist_unit = "mi", st.bottom = FALSE, st.color = "white", transform = TRUE, model = "WGS84") +
  theme(legend.position="none",
        axis.text = element_text(size=11)) +
  labs(x="", y="")

ggsave("~/github/mysosoup/figures/Figure1_rawMap_BothSites.png")
ggsave("~/github/mysosoup/figures/Figure1_rawMap_BothSites.svg",
       height = 15, width = 15, units = "cm")

############ 

hick_sat <- get_googlemap(center=c(lon = -89.072, lat = 37.339), zoom=15, maptype = "satellite")
ggmap(hick_sat) +
  geom_point(aes(x=lon, y=lat, color=Area), data=sitedata %>% filter(Area=="Hickory"), size=2.5) +
  scale_color_manual(values=c("pink")) +
  scalebar(x.min = -89.06, x.max = -89.072, y.min = 37.33, y.max = 37.332,
           dist = 0.5, dist_unit = "mi", st.bottom = FALSE, st.color = "white", transform = TRUE, model = "WGS84") +
  theme(legend.position="none",
        axis.text = element_text(size=11)) +
  labs(x="", y="")

## save as "Figure1_rawMap_HickOnly"; export at 750x750

############ 

egnr_sat <- get_googlemap(center=c(lon = -89.095, lat = 37.295), zoom=15, maptype = "satellite")
ggmap(egnr_sat) +
  geom_point(aes(x=lon, y=lat, color=Area), data=sitedata %>% filter(Area=="Egner"), size=2.5) +
  scale_color_manual(values=c("gold2")) +
  scalebar(x.min = -89.083, x.max = -89.12, y.min = 37.286, y.max = 37.287,
           dist = 0.5, dist_unit = "mi", st.bottom = FALSE, st.color = "white", transform = TRUE, model = "WGS84") +
  theme(legend.position="none",
        axis.text = element_text(size=11)) +
  labs(x="", y="")

## save as "Figure1_rawMap_EgnerOnly"; export at 750x750


###################################################################
## unused code
###################################################################


## read in sample metadata to import the number of samples collected at each site, so we don't have to add a table about monthly collections?

df <- read_csv("https://github.com/devonorourke/tidybug/raw/master/data/text_tables/all.filtmethods.df.csv.gz") %>% 
  filter(SampleType == "sample")

sitedata <- metadata %>% 
  filter(SampleID %in% df$SeqID) %>% 
  group_by(Site, SiteLat, SiteLong) %>% 
  tally()

sitedata$SiteLat <- ifelse(sitedata$Site == "FAR", 44.5884, sitedata$SiteLat)
sitedata$SiteLong <- ifelse(sitedata$Site == "FAR", -69.5986, sitedata$SiteLong)

sitedata <- sitedata %>% 
  rename(lon = SiteLong, lat = SiteLat) %>% 
  mutate(PlotLabel = paste(Site, n, sep = "\n"))

## import map layers from Google
basemap <- get_stamenmap(bbox = c(left = -72.5, bottom = 42.5, right = -69.5, top = 45),
                         zoom = 9,maptype = "toner-background")

ggmap(basemap) +
  geom_point(aes(x=lon, y=lat), color="firebrick", data=sitedata, size=2.5) +
  geom_label_repel(aes(x=lon, y=lat, label=PlotLabel), data=sitedata, size=3)


## save as "FigureS1_collectionMap"; export at 650x650

