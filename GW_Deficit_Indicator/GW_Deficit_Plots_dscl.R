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
GGDI.out = fread(paste0(pathOut, 'GGDI_output_dscl_230128.csv'))

ggdi.r = rasterFromXYZ(GGDI.out[,.(lon, lat, Def.19_20_neg15)])

#####Other Rasters/Data
#WB regions raster
globe.r = 
  fasterize(aq.v, fishnet.r)

#Downscaled regions
dscl = 
  rasterFromXYZ(GGDI.out[,.(lon, lat, Def.19_20_neg15)]) %>%
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



var2plot = c('Def.19_20_neg15', 'Def.19_20_neg1', 'Def.bin24_150',
             'Def.bin24_100', 'neg_sig')


i = 1

for(i in 1:length(var2plot)){
  
  cur_var = var2plot[i]
  plot_var = str_replace(cur_var, '\\.', '_')
  
  if(cur_var %in% c('Def.19_20_neg15', 'Def.19_20_neg1')) {
    label = 'Deficit 19/20\n(downscaled)'
  } else if (cur_var %in% c('Def.bin24_150', 'Def.bin24_100')) {
    label = 'Deficit 02/20\n(downscaled)'
  } else if (cur_var == 'neg_sig') {
    label = 'Negative Trend\n (p<0.05)'
  }
  
  cur_rast = 
    GGDI.out %>% 
    dplyr::select('lon', 'lat', all_of(cur_var)) %>%
    rasterFromXYZ()
  
  crs(cur_rast) = crs(wb_regions)
  
  cur_rast[values(cur_rast)==0,] = NA
  
  cur_rast = 
    crop(cur_rast, extent(-180, 180, -60, 60)) %>%
    as.data.frame(xy=T, na.rm = T) %>%
    mutate(deficit=factor(all_of(cur_var))) %>%
    as.data.table()
  
  levels(cur_rast$deficit) = label
  
  ## SETUP FOR BAR AND BOXPLOTS
  scale.aq <- scale_fill_manual(values=c("#44546a", "#70ad47", "#b7ff4b", '#ffc000'), name = "Aquifer type")
  theme.blank <-   theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank())
  
  plot1 = 
    aq.v %>%
    drop_na() %>%
    ggplot() +
    geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
    scale.aq +
    geom_sf(data = problem, fill = '#ffc000', alpha = 0.1, lwd = 0, show.legend = FALSE) + 
    new_scale_fill() +
    geom_tile(data = cur_rast, aes(x = x, y = y, fill = deficit), alpha = 0.85) +
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
  
  
  ggsave(paste0(pathOut, 'Figures/', plot_var, '_aqtyp_GB', '.png'), plot=plot1,
         scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')
  
  
  plot2 = 
    aq.v.dscl %>%
    drop_na() %>%
    ggplot() +
    geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
    scale.aq +
    new_scale_fill() +
    geom_tile(data = cur_rast, aes(x = x, y = y, fill = deficit), alpha = 0.85) +
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
  
  ggsave(paste0(pathOut, 'Figures/', plot_var, '_aqtyp_SR', '.png'), plot=plot2,
         scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')
  
  
  
  plot3 = 
    ggplot() +
    geom_sf(data = studyRegion, fill = NA, size = 0.8) +
    new_scale_fill() +
    geom_tile(data = cur_rast, aes(x = x, y = y, fill = deficit), alpha = 0.8) +
    scale_fill_manual(values='#d6604d', name = 'Deficit') +
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
  
  ggsave(paste0(pathOut, 'Figures/', plot_var, '_TBA', '.png'), plot=plot3,
         scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')
  
}




