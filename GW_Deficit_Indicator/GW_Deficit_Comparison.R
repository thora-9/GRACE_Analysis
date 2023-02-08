#Compare GW deficit indicators - downscaled vs non-downscaled

library(tidync)
library(data.table)
library(tidyverse)
library(lubridate)
library(zyp)
library(sf)
library(raster)
library(terra)
library(RColorBrewer)
library(rasterVis)
library(xts)
library(stars)
library(fasterize)
library(foreign)
library(RcppRoll)
library(tmap)

############################################################################################
###Load  the paths
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"
pathOut =  "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"
pathIn = "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"
pathData = '/Users/tejasvi/Dropbox/gwflagship_GRACEdownscaling/Downscaled TWS_GWS v2/'
pathTyp = '/Users/tejasvi/Dropbox/WB/Typology/'


mymerge = function(x,y) {merge.data.table(x,y,all=F,
                                          by.x = c('lat', 'lon'),
                                          by.y = c('Lat', 'Lon'))}

############################################################################################

#Load the studyregions
med.v = st_read(paste0(pathData, 'studyregions/med.shp'))
sa.v = st_read(paste0(pathData, 'studyregions/sa.shp'))
saf.v = st_read(paste0(pathData, 'studyregions/saf.shp'))

studyRegion = 
  st_union(med.v, sa.v) %>%
  st_union(saf.v)

####Load World Regions
wb_regions = 
  st_read(paste0(proj_dir, "Spatial Files/WB_Regions/WB_countries_Admin0_10m.shp")) %>%
  dplyr::select(WB_NAME, ISO_A2, ISO_A3, ISO_N3, TYPE, REGION_WB) %>%
  filter(TYPE != 'Dependency') %>%
  st_make_valid()

wb_region_dissolved = 
  wb_regions %>%
  group_by(REGION_WB) %>% 
  summarize(geometry = st_union(geometry)) %>%
  mutate(ID = row_number()) %>%
  mutate(REGION_WB = as.character(REGION_WB)) %>%
  #To combine MENA and SSA into one (to avoid dangling cross-region TBA)
  mutate(REGION_WB = ifelse(REGION_WB %in% c('Middle East & North Africa','Sub-Saharan Africa'),
                            'Middle East & Africa', REGION_WB)) %>%
  group_by(REGION_WB) %>% 
  summarize(geometry = st_union(geometry)) 

#######
#Aquifer Typology
## Vector version
aq.v <-
  st_read(paste0(pathTyp, "aqtyp_dissolved.gpkg")) %>%
  # dplyr::select(aqtyp) %>%
  # st_union()
  mutate(aqtyp=factor(aqtyp, levels=c("Major Alluvial","Complex","Karstic","Local/Shallow")))

problem =
  aq.v[84800,]


## transboundary aquifers
trans.sub = c('AF058','AF064','AF056', 'AF001','AF063', #'AF054', Remove Volta
              'AS126','AS150','AS079','AS080','AS089','AS090','AS091',
              'AS150','S021','S015','C007')

trans <-
  st_read(paste0(pathTyp, 'IGRAC_Transboundary/2021_TBA_GGIS_utf8_VALID.shp')) %>%
  filter(CODE_2021 %in% trans.sub) %>%
  st_make_valid() %>%
  mutate(area = st_area(.) %>% as.numeric()) %>%
  mutate(area_ha = area*0.0001) %>%
  st_join(wb_region_dissolved[,'REGION_WB'], left = T) %>%
  filter(REGION_WB %in% c('South Asia', 'Middle East & Africa')) %>%
  distinct()

trans.all <- 
  st_read(paste0(pathTyp, 'IGRAC_Transboundary/2021_TBA_GGIS_utf8_VALID.shp')) %>%
  st_make_valid() %>%
  st_join(wb_region_dissolved[,'REGION_WB'], left = T) %>%
  filter(REGION_WB %in% c('South Asia', 'Middle East & Africa')) %>%
  dplyr::select(-REGION_WB) %>%
  distinct() 

#######
files = list.files(pathIn, '230128', full.names = T)

files.names = 
  list.files(pathIn, '230128', full.names = F) 
files.names = gsub("\\..*","",files.names)

#Load the datasets
ggdi.data = fread(files)

#Merge the datasets
ggdi.data = 
  ggdi.data %>%
  mutate(sumDef = ifelse(Def.bin24_100  == 0 & neg_sig == 0, 1, 
                         ifelse(Def.bin24_100  == 1 & neg_sig == 0, 2,
                                ifelse(Def.bin24_100  == 0 & neg_sig == 1, 3, 4))))



############################################################################################

#Overall
table(ggdi.data$Def.bin24_100,
      ggdi.data$neg_sig,
      dnn = c("Def_10_20", "neg_sig"), useNA = 'no')

#Aquifer Type
tab.aqtyp = 
  table(ggdi.data$Def.bin24_100,
        ggdi.data$aqtyp_max,
        dnn = c("Def_10_20", "aquifer type"), useNA = 'no') %>% as.data.frame() %>%
  group_by(aquifer.type) %>%
  mutate(total = sum(Freq)) %>%
  mutate(prop = Freq/total)

temp.rast = 
  ggdi.data %>% 
  dplyr::select('lon', 'lat', 'neg_sig') %>%
  rasterFromXYZ()

#Percent negative trends in each TBA
TBA_neg = 
  raster::extract(temp.rast, 
                  trans.all)

len1 = sapply(TBA_neg, function(x){length(x)})
sumNA = sapply(TBA_neg, function(x){sum(is.na(x))})
sumNA.n = sapply(TBA_neg, function(x){sum(!is.na(x))})
sum = sapply(TBA_neg, function(x){sum(x, na.rm = T)})


TBA_neg = 
  cbind(trans.all, len1, sumNA, sumNA.n, sum) %>%
  mutate(out.SR = sumNA/len1) %>%
  filter(out.SR != 1) %>%
  mutate(atleast1 = ifelse(sum>0, 1, 0))

100*sum(TBA_neg$sum>0)/nrow(TBA_neg)
#st_write(TBA_neg, paste0(pathOut, 'TBA_neg.shp'), append = F)

cur_rast = 
  ggdi.data %>% 
  dplyr::select('lon', 'lat', 'sumDef') %>%
  rasterFromXYZ()

crs(cur_rast) = crs(wb_regions)

cur_rast[values(cur_rast)==1,] = NA

cur_rast = 
  crop(cur_rast, extent(-180, 180, -60, 60)) %>%
  as.data.frame(xy=T, na.rm = T) %>%
  mutate(deficit=factor(sumDef)) %>%
  as.data.table()

levels(cur_rast$deficit) = c("Deficit 10/20", "Negative Trend", "Trend & Deficit")


#Clip the aquifer typlogy raster
aq.v.clipped = 
  aq.v[-84800, ] %>%
  st_make_valid()

aq.v.dscl = 
  st_crop(aq.v.clipped, temp.rast)


cls <- data.frame(id=1:4, cover=c("No Deficit", 
                                  "Deficit 10/20 only",
                                  "Deficit 19/20 only",
                                  "Deficit - Both"))

############################################################################################
## SETUP FOR BAR AND BOXPLOTS
scale.aq <- scale_fill_manual(values=c("#44546a", "#70ad47", "#b7ff4b", '#ffc000'), name = "Aquifer type")
theme.blank <-   theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())


plot1 = 
  aq.v.dscl %>%
  drop_na() %>%
  ggplot() +
  geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
  scale.aq +
  new_scale_fill() +
  geom_tile(data = cur_rast, aes(x = x, y = y, fill = deficit), alpha = 0.85) +
  scale_fill_manual(values=c('#787878', '#f7918d', '#d6604d'), name = 'Deficit') +
  new_scale_fill() +
  geom_sf(data = studyRegion, fill = NA, size = 0.8) +
  theme_bw() + ylim(-34.83427, 37.5) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathOut, 'Figures/Def_comparison_aqtyp_230204', '.png'), plot=plot1,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')

plot2 = 
  ggplot() +
  geom_sf(data = studyRegion, fill = NA, size = 0.8) +
  new_scale_fill() +
  geom_tile(data = cur_rast, aes(x = x, y = y, fill = deficit), alpha = 0.8) +
  scale_fill_manual(values=c('#787878', '#f7918d', '#d6604d'), name = 'Deficit') +
  theme_bw() + ylim(-34.83427, 37.5) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathOut, 'Figures/Def_comparison_230204', '.png'), plot=plot2,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')


plot3 = 
  ggplot() +
  geom_sf(data = studyRegion, fill = NA, size = 0.8) +
  new_scale_fill() +
  geom_tile(data = cur_rast, aes(x = x, y = y, fill = deficit), alpha = 0.8) +
  scale_fill_manual(values=c('#787878', '#f7918d', '#d6604d'), name = 'Deficit') +
  new_scale_fill() +
  geom_sf_pattern(data = trans, 
                  aes(pattern = 'Transboundary\n Aquifer'), colour = 'firebrick4',
                  fill = NA, pattern_colour  = 'black', pattern_alpha = 0.7, 
                  pattern_density = 0.6, pattern_fill = NA,
                  pattern_spacing = 0.005) +
  geom_label_repel(data = trans,
                   aes(label = CODE_2021, geometry = geometry),
                   stat = "sf_coordinates",
                   min.segment.length = 0) +
  scale_pattern_manual(name = NA, values = 'circle') +
  theme_bw() + ylim(-34.83427, 37.5) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathOut, 'Figures/Def_comparison_TBA_230204', '.png'), plot=plot3,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')


plot4 = 
  ggplot() +
  geom_sf(data = studyRegion, fill = NA, size = 0.8) +
  new_scale_fill() +
  geom_tile(data = cur_rast, aes(x = x, y = y, fill = deficit), alpha = 0.8) +
  scale_fill_manual(values=c('#787878', '#f7918d', '#d6604d'), name = 'Deficit') +
  new_scale_fill() +
  geom_sf_pattern(data = trans.all, 
                  aes(pattern = 'Transboundary\n Aquifer'), colour = 'firebrick4',
                  fill = NA, pattern_colour  = 'black', pattern_alpha = 0.7, 
                  pattern_density = 0.6, pattern_fill = NA,
                  pattern_spacing = 0.005) +
  geom_label_repel(data = trans.all,
                   aes(label = CODE_2021, geometry = geometry),
                   stat = "sf_coordinates",
                   min.segment.length = 0) +
  scale_pattern_manual(name = NA, values = 'circle') +
  theme_bw() + ylim(-34.83427, 37.5) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathOut, 'Figures/Def_comparison_TBAall_230204', '.png'), plot=plot4,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')
