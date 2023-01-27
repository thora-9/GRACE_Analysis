#Produce groundwater deficit indicator datasets (downscaled)

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
library(haven)
library(foreign)
library(stars)
library(ggnewscale)
library(ggpattern)
library(gridExtra)
library(ggrepel)

############################################################################################
pathTyp = '/Users/tejasvi/Dropbox/WB/Typology/'
pathIn = '/Users/tejasvi/Dropbox/gwflagship_typologies/'
proj_dir = "/Users/tejasvi/Dropbox/WB/GRACE_Ensemble/"
pathOut =  "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"
pathData = '/Users/tejasvi/Dropbox/gwflagship_GRACEdownscaling/Downscaled TWS_GWS v2/'

############################################################################################

####Load World Regions
wb_regions = 
  st_read(paste0(proj_dir, "Spatial Files/WB_Regions/WB_countries_Admin0_10m.shp")) %>%
  dplyr::select(WB_NAME, ISO_A2, ISO_A3, ISO_N3, TYPE, REGION_WB) %>%
  filter(TYPE != 'Dependency') %>%
  st_make_valid()

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

#######
#Use the fishnet to expand the extent of the Downscaled ouputs
#Load the fishnet
fishnet = 
  st_read('/Users/tejasvi/Dropbox/WB/Fishnet_halfdegree/global_fishnet.shp')
fishnet.r = 
  fishnet %>%
  st_drop_geometry() %>%
  as.data.table()
fishnet.r = 
  rasterFromXYZ(fishnet.r[,.(Lon, Lat)])

crs(fishnet.r) = crs(wb_regions)

values(fishnet.r) = 1

#######
#Groundwater Deficit
GGDI.out = fread(paste0(pathOut, 'GGDI_output_dscl_230103.csv'))

ggdi.r = rasterFromXYZ(GGDI.out[,.(lon, lat, Def.19_20)])

#######
#2019-2020 Indicator
ggdi.19 = 
  rasterFromXYZ(GGDI.out[,.(lon, lat, Def.19_20)])

crs(ggdi.19) = crs(wb_regions)

ggdi.19[values(ggdi.19)==0,] = NA

ggdi.19 = 
  crop(ggdi.19, extent(-180, 180, -60, 60)) %>%
  as.data.frame(xy=T, na.rm = T) %>%
  mutate(deficit=factor(Def.19_20)) %>%
  as.data.table()
levels(ggdi.19$deficit) = 'Deficit 19/20\n(downscaled)'

#######
#2002-2020 Indicator
ggdi.02 = 
  rasterFromXYZ(GGDI.out[,.(lon, lat, Def.bin24_150)])

crs(ggdi.02) = crs(wb_regions)

ggdi.02[values(ggdi.02)==0,] = NA

ggdi.02 = 
  crop(ggdi.02, extent(-180, 180, -60, 60)) %>%
  as.data.frame(xy=T, na.rm = T) %>%
  mutate(deficit=factor(Def.bin24_150)) %>%
  as.data.table()
levels(ggdi.02$deficit) = 'Deficit 02/20\n(downscaled)'

#WB regions raster
globe.r = 
  fasterize(aq.v, fishnet.r)

#Downscaled regions
dscl = 
  rasterFromXYZ(GGDI.out[,.(lon, lat, Def.19_20)]) %>%
  extend(fishnet.r)
crs(dscl) = crs(wb_regions)
dscl[!is.na(values(dscl)), ] = 1

aq.v.clipped = 
  aq.v[-84800, ] %>%
  st_make_valid()

aq.v.dscl = 
  st_crop(aq.v.clipped, ggdi.r)

#Non-downscaled regions
non.dscl = 
  globe.r
values(non.dscl) = NA
non.dscl[values(globe.r)==1,] = 1
non.dscl[values(dscl)==1,] = NA

non.dscl.df = 
  non.dscl %>%
  as.data.frame(xy=T, na.rm = T) 

non.dscl.dscl = 
  crop(non.dscl, ggdi.r) %>%
  as.data.frame(xy=T, na.rm = T) 


#Load the studyregions
med.v = st_read(paste0(pathData, 'studyregions/med.shp'))
sa.v = st_read(paste0(pathData, 'studyregions/sa.shp'))
saf.v = st_read(paste0(pathData, 'studyregions/saf.shp'))

studyRegion = 
  st_union(med.v, sa.v) %>%
  st_union(saf.v)
############################################################################################
#Plots
############################################################################################
## SETUP FOR BAR AND BOXPLOTS
scale.aq <- scale_fill_manual(values=c("#44546a", "#70ad47", "#b7ff4b", '#ffc000'), name = "Aquifer type")
theme.blank <-   theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())


#Global Map

plot2 = 
  aq.v %>%
  drop_na() %>%
  ggplot() +
  geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
  scale.aq +
  geom_sf(data = problem, fill = '#ffc000', alpha = 0.1, lwd = 0, show.legend = FALSE) + 
  new_scale_fill() +
  geom_tile(data = ggdi.19, aes(x = x, y = y, fill = deficit), alpha = 0.85) +
  scale_fill_manual(values='#d6604d', name = 'Deficit') +
  new_scale_fill() +
  geom_tile(data = non.dscl.df, aes(x = x, y = y), fill = 'grey', alpha = 0.95) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

# ggsave(paste0(pathOut, 'Figures/GW_Deficit_19_20_aqtyp_dscl', '.png'), plot=plot2,
#        scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')



plot3 = 
  aq.v %>%
  drop_na() %>%
  ggplot() +
  geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
  scale.aq +
  geom_sf(data = problem, fill = '#ffc000', alpha = 0.1, lwd = 0, show.legend = FALSE) + 
  new_scale_fill() +
  geom_tile(data = ggdi.02, aes(x = x, y = y, fill = deficit), alpha = 0.8) +
  scale_fill_manual(values='#d6604d', name = 'Deficit') +
  new_scale_fill() +
  geom_tile(data = non.dscl.df, aes(x = x, y = y), fill = 'grey', alpha = 0.95) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

# ggsave(paste0(pathOut, 'Figures/GW_Deficit_02_20_aqtyp_dscl', '.png'), plot=plot3,
#        scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')


#########################
plot2 = 
  aq.v.dscl %>%
  drop_na() %>%
  ggplot() +
  geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
  scale.aq +
  new_scale_fill() +
  geom_tile(data = ggdi.19, aes(x = x, y = y, fill = deficit), alpha = 0.85) +
  scale_fill_manual(values='#d6604d', name = 'Deficit') +
  new_scale_fill() +
  geom_sf(data = studyRegion, fill = NA, size = 0.8) +
  #geom_tile(data = non.dscl.dscl, aes(x = x, y = y), fill = 'white', alpha = 0.95) +
  theme_bw() + ylim(-34.83427, 37.5) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathOut, 'Figures/GW_Deficit_19_20_aqtyp_dscl_ALT', '.png'), plot=plot2,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')



plot3 = 
  aq.v.dscl %>%
  drop_na() %>%
  ggplot() +
  geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
  scale.aq +
  new_scale_fill() +
  geom_tile(data = ggdi.02, aes(x = x, y = y, fill = deficit), alpha = 0.8) +
  scale_fill_manual(values='#d6604d', name = 'Deficit') +
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

ggsave(paste0(pathOut, 'Figures/GW_Deficit_02_20_aqtyp_dscl_ALT', '.png'), plot=plot3,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')
