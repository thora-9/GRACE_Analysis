#Produce groundwater deficit indicator datasets

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
###Load all the GWS datasets
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"
pathOut =  "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"
pathIn = '/Users/tejasvi/Dropbox/gwflagship_typologies/'

############################################################################################
####Load World Regions
wb_regions = 
  st_read(paste0(proj_dir, "Spatial Files/WB_Regions/WB_countries_Admin0_10m.shp")) %>%
  dplyr::select(WB_NAME, ISO_A2, ISO_A3, ISO_N3, TYPE, REGION_WB) %>%
  filter(TYPE != 'Dependency') %>%
  st_make_valid()


#load the non-downscaled GRACE data
gws.OG = 
  fread(paste0(proj_dir, 
               'Outputs/GWS/half-degree/GRACE_GWS_3Mascons_2002_2020_BSL2017_05degree.csv')) %>%
  dplyr::rename(cell_id_OG = cell_id) %>%
  #st_drop_geometry() %>%
  mutate(year = substr(ym, 1, 4),
         month = substr(ym, 6, 7),
         cell_id = paste0(Lon, Lat))

gws.unique = 
  gws.OG %>% 
  dplyr::select(-year,-month, -Shape_Leng, -Shape_Area, -GWS_ensemble, -ym) %>%
  distinct() %>%
  rownames_to_column() %>%
  as.data.table()

  
plot1 = 
  rasterFromXYZ(gws.unique[,.(Lon, Lat, ISO_N3)])

#Load the fishnet
fishnet = 
  st_read('/Users/tejasvi/Dropbox/WB/Fishnet_halfdegree/global_fishnet.shp')

# fishnet.centroid = 
#   fishnet %>%
#   filter(!row_number() %in% c(18901, 23763, 26074, 42334, 31737)) %>% #Invalid polygon
#   st_make_valid() %>%
#   st_centroid()

fishnet.r = 
  fishnet %>%
  st_drop_geometry() %>%
  as.data.table()

fishnet.r = 
  rasterFromXYZ(fishnet.r[,.(Lon, Lat)])

crs(fishnet.r) = crs(wb_regions)

#Load aquifer typology data
typ <- 
  read.dta(file=paste0(pathIn, "data_outputs/aqtyp_gwresource_grid05deg.dta"), convert.factors = TRUE) %>% 
  as_tibble %>%
  mutate(aqtyp=factor(aqtyp_max, levels=c("Major Alluvial","Complex","Karstic","Local/Shallow"))) %>%
  as.data.table()

############################################################################################
#GGDI steps
############################################################################################
GGDI = 
  gws.OG %>%
  filter(year<2021) %>%
  merge(gws.unique[,.(Lat, Lon, rowname)], by = c('Lat', 'Lon'), all.x=T) %>%
  rename(GWS = GWS_ensemble) %>%
  group_by(rowname, month) %>% #Group by month and cell to get climatology
  mutate(GWS.month.mean = mean(GWS, na.rm = T)) %>% 
  rowwise() %>%
  mutate(GWS.climatology = GWS - GWS.month.mean) %>%
  group_by(rowname) %>%
  mutate(grid.mean = mean(GWS.climatology, na.rm = T),
         grid.sd = sd(GWS.climatology, na.rm = T)) %>%
  mutate(GWS.deficit = (GWS.climatology-grid.mean)/grid.sd) %>%
  mutate(GWS.def.roll12 = roll_mean(GWS.deficit, 12, align = 'right', fill = NA),
         GWS.def.roll18 = roll_mean(GWS.deficit, 18, align = 'right', fill = NA),
         GWS.def.roll24 = roll_mean(GWS.deficit, 24, align = 'right', fill = NA)) %>%
  as.data.table()
  
View(GGDI[cell_id == '25.2568.25'])


#Basic Indicator -  Is the GW Deficit value between 2019-2020 less than -1.5
GGDI.binary1 = 
  GGDI %>%
  filter(year>=2019) %>%
  dplyr::select(rowname, GWS.deficit) %>%
  group_by(rowname) %>%
  summarise(GWS.deficit.mean = mean(GWS.deficit, na.rm = T)) %>% 
  mutate(Def.19_20 = ifelse(GWS.deficit.mean > -1.5, 0, 1)) %>%
  as.data.table()

#Basic Indicator -  Is the GW Deficit value between 2013-2014 less than -1.5 (to coincide with Esha's dataset)
GGDI.binary1b = 
  GGDI %>%
  filter(year>=2013 & year<2015) %>%
  dplyr::select(rowname, GWS.deficit) %>%
  group_by(rowname) %>%
  summarise(GWS.deficit.mean = mean(GWS.deficit, na.rm = T)) %>% 
  mutate(Def.13_14 = ifelse(GWS.deficit.mean > -1, 0, 1)) %>%
  as.data.table()

#Indicator2 -  Is there any period between 2002 and 2017 where the GW Deficit value is less than -1.5
GGDI.binary2 = 
  GGDI %>%
  filter(year<=2017) %>%
  dplyr::select(rowname, GWS.def.roll12, GWS.def.roll18, GWS.def.roll24) %>%
  mutate(GWS.binary12 = ifelse(GWS.def.roll12 > -1.25, 0, 1),
         GWS.binary18 = ifelse(GWS.def.roll18 > -1.5, 0, 1),
         GWS.binary24_150 = ifelse(GWS.def.roll24 > -1.5, 0, 1),
         GWS.binary24_125 = ifelse(GWS.def.roll24 > -1.25, 0, 1)) %>%
  group_by(rowname) %>%
  summarise(Def.total12 = sum(GWS.binary12, na.rm = T),
            Def.total18 = sum(GWS.binary18, na.rm = T),
            Def.total24_150 = sum(GWS.binary24_150, na.rm = T),
            Def.total24_125 = sum(GWS.binary24_125, na.rm = T)) %>%
  mutate(Def.bin12 = ifelse(Def.total12>0, 1, 0),
         Def.bin18 = ifelse(Def.total18>0, 1, 0),
         Def.bin24_150 = ifelse(Def.total24_150>0, 1, 0),
         Def.bin24_125 = ifelse(Def.total24_125>0, 1, 0)) %>%
  as.data.table() 



#Merge indicators
GGDI.out = 
  GGDI.binary1 %>%
  merge(GGDI.binary1b[,.(rowname, Def.13_14)], by = 'rowname', all.x = T) %>%
  merge(GGDI.binary2, by = 'rowname', all.x = T) %>%
  merge(gws.unique, by = 'rowname', all.x = T) %>%
  merge(typ[,c('lat', 'lon', 'aqtyp_max')], by.x = c('Lat', 'Lon'), by.y = c('lat', 'lon')) %>%
  #Basically, making sure that hotspots in 2019 remain hotspots in the rolling mean indicator
  mutate(Def.bin24_150 = ifelse(Def.19_20==1, 1, Def.bin24_150)) 
  
GGDI.out.sub = 
  GGDI.out %>%
  dplyr::select(1:6, Def.total24_150, Def.bin24_150, Id:aqtyp_max)


############################################################################################
#Plots
############################################################################################

#####
plot1 = 
  rasterFromXYZ(GGDI.out[,.(Lon, Lat, Def.19_20)])

crs(plot1) = crs(fishnet.r)

plot1 = 
  crop(plot1, extent(-180, 180, -60, 60)) %>%
  rast()

cls <- data.frame(id=0:1, cover=c("No Deficit", "Deficit"))
levels(plot1) <- cls


#color_pal = c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')

#Save the plot for future reference
plot_name =
  paste0("~/Dropbox/WB/GRACE-Deficit/Figures/", "GW_Deficit_19_20.png")

png(plot_name, width = 1250, height = 500)

plot(plot1,
     col = c('#4393c3', '#d6604d'), plg=list(cex=1.2))

dev.off()

# 
# plot2 = 
#   rasterFromXYZ(GGDI.out[,.(Lon, Lat, Def.total)])
# 
# crs(plot2) = crs(fishnet.r)
# 
# plot2 = 
#   crop(plot2, extent(-180, 180, -60, 60))
# 
# 
# color_pal = c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')
# 

# plot(plot2,
#      legend.args = list(text = 'Depletion', side = 4, 
#                         font = 2, line = 2.5, cex = 0.8)) 
# 

#####

plot3 = 
  rasterFromXYZ(GGDI.out[,.(Lon, Lat, Def.bin24_150)])

crs(plot3) = crs(fishnet.r)

plot3 = 
  crop(plot3, extent(-180, 180, -60, 60)) %>%
  rast()

cls <- data.frame(id=0:1, cover=c("No Deficit", "Deficit"))
levels(plot3) <- cls

color_pal = c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')


#Save the plot for future reference
plot_name =
  paste0("~/Dropbox/WB/GRACE-Deficit/Figures/", "GW_Deficit_24_month_binary.png")

png(plot_name, width = 1250, height = 500)

plot(plot3,
     col = c('#4393c3', '#d6604d'), plg=list(cex=1.2))

dev.off()



############################################################################################
#Distribution of hotspots by aq typ
table(GGDI.out[!aqtyp_max %in% c("", NA) & Def.19_20==1,]$aqtyp_max) %>% prop.table()
table(GGDI.out[!aqtyp_max %in% c("", NA) & Def.bin24_150==1,]$aqtyp_max, useNA = 'no') %>% prop.table()

#Overall
table(GGDI.out$Def.19_20,
      GGDI.out$Def.bin24_150, 
      dnn = c("Def.19_20", "Def.bin24_150"), useNA = 'no') 

#Without NA
table(GGDI.out[!WB_REGION %in% c("Other"),]$Def.19_20,
      GGDI.out[!WB_REGION %in% c("Other"),]$Def.bin24_150, 
      dnn = c("Def.19_20", "Def.bin24_150"), useNA = 'no') %>% prop.table()


############################################################################################
#Write Files
############################################################################################


pathOut =  "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"

fwrite(GGDI.out.sub, paste0(pathOut, 'GGDI_output_nonDownscaled_230103.csv'))

# 
st_write(GGDI.fishnet,
         paste0(pathOut, 'GWS_Deficit_nonDownscaled_05degree_230103.shp'),
         delete_layer = T)

