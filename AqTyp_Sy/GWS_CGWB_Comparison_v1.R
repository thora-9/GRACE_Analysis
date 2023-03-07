#Link aquifer types in India with specific yield from CGWB

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
library(ggnewscale)
library(readr)
library(strex)

############################################################################################
###Load all the GWS datasets
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"
pathOut =  "/Users/tejasvi/Dropbox/WB/Downscaling GRACE/"
pathIn = "/Users/tejasvi/Dropbox/WB/India_Aquifer_Sy/"
pathData = '/Users/tejasvi/Dropbox/gwflagship_GRACEdownscaling/Downscaled GWS v3/'

############################################################################################

gws_mean_05 = 
  tidync(paste0(pathData, 'GWS_WB_025deg.nc')) %>%
  hyper_tibble() %>% as.data.table()

#Create a date sequence
date.seq = 
  seq(as.Date("2003/2/1"), as.Date("2021/9/1"), "month") %>%
  as.data.table() %>%
  rownames_to_column() %>%
  .[, rowname := as.integer(rowname)]
colnames(date.seq) = c('rowname', 'date')

#Merge the dates using the rowname column
gws_mean_05 = 
  gws_mean_05 %>%
  merge(date.seq, by.x = 'time', by.y = 'rowname', all.x = T) %>%
  .[, ':='(yearmon = substr(date, 1, 7))] %>%
  dplyr::rename(GWSA_mean = GWS_WB_025deg)

#Convert to wide
gws_05_wide = 
  gws_mean_05 %>% 
  dplyr::select(GWSA_mean, lat, lon, date) %>%
  tidyr::spread(key = date, value = GWSA_mean) %>%
  mutate(cell_id = paste0(lon, lat)) %>%
  as.data.table()

#Convert to spatial object
gws_05.v = 
  st_as_sf(gws_05_wide, coords = c("lon", "lat"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs")

####Load World Regions
wb_regions = 
  st_read(paste0(proj_dir, "Spatial Files/WB_Regions/WB_countries_Admin0_10m.shp")) %>%
  dplyr::select(WB_NAME, ISO_A2, ISO_A3, ISO_N3, TYPE, REGION_WB) %>%
  filter(TYPE != 'Dependency') %>%
  st_make_valid()

#Merge country data with GRACE
gws_05.v = 
  gws_05.v %>%
  st_make_valid() %>%
  st_join(wb_regions)

#Need to isolate south asia for plots
gws_05_sa = 
  gws_05.v %>% filter(as.character(WB_NAME) == 'India') 

gws_temp = gws_05_wide %>% filter(cell_id %in% gws_05_sa$cell_id)

gws.raster = 
  rasterFromXYZ(gws_temp[,.(lon, lat, `2003-02-01`)])
#values(gws.raster) = 1


############################################################################################
#Need to essentially use the gws.raster to create a fishnet polygon
#which will be used to compare the output from GRACE and in-situ
rasterTemplate = gws.raster
vectorTemplate = 
  rasterTemplate %>%
  st_as_stars() %>%
  st_as_sf(as_points = FALSE, crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  rownames_to_column()

st_crs(vectorTemplate) = "+proj=longlat +datum=WGS84 +no_defs"

############################################################################################

#Load the in-situ CGWB data
cgwb.dug = 
  fread(paste0(pathIn, 'CGWB_w_specificYield.csv')) %>%
  dplyr::select(-(`May 1996`:`Nov 2002`)) %>%
  mutate(Nas = rowSums(is.na(.[,6:62]))) %>% #Keep only the rows after 2003
  dplyr::filter(Nas<(56-24) & SITE_TYPE == 'Dug Well') %>% #Filter out wells with less than 6 years of data
  st_as_sf(coords = c("LON", "LAT"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_join(vectorTemplate[,c('rowname')]) %>%
  dplyr::select(-V1, -STATE, -DISTRICT, -SITE_TYPE, -Nas)

colname.ref = colnames(cgwb.dug)

cgwb.long = 
  cgwb.dug %>%
  st_drop_geometry() %>% as.data.table() %>%
  melt(id.vars    = c("WLCODE", "rowname", "sy_final"),
       measure.vars  = colname.ref[2:58],
       variable.name = "yearmon",
       value.name    = "measured") %>%
  group_by(WLCODE) %>%
  mutate(mean_all = mean(measured, na.rm = T)) %>% 
  rowwise() %>%
  mutate(gwl.anomaly = -(measured - mean_all)) %>% #Subtract the all-time mean to obtain a anomaly; take negative to get depth below surface
  mutate(gws.anomaly = gwl.anomaly*(sy_final/100))


#Now summarize by grid cell
cgwb.grid = 
  cgwb.long %>% 
  group_by(rowname, yearmon) %>%
  dplyr::summarise(measured.mean = mean(gws.anomaly, na.rm = T)) %>%
  as.data.table() %>%
  mutate(yearmon = as.character(yearmon))


#Link grid-cell id to grace data
gws_05.SA = 
  gws_05.v %>% 
  filter(as.character(WB_NAME) == 'India') %>%
  st_join(vectorTemplate[,c('rowname')]) %>%
  st_drop_geometry() %>% as.data.table() %>%
  data.table::melt(id.vars = c("cell_id", "rowname"),
                   measure.vars = 3:212,
                   variable.name = "ymd", value.name = "GWSA_mean") %>%
  mutate(yearmon = as.yearmon(as.Date(ymd)))

#Create a merged dataset
gws.cgwb.merged = 
  gws_05.SA %>% dplyr::select(rowname, yearmon, ymd, GWSA_mean) %>%
  mutate(yearmon = as.character(yearmon),
         ymd = as.Date(ymd)) %>%
  merge(cgwb.grid, by = c('yearmon', 'rowname')) %>%
  arrange(rowname, ymd)

#Assess the correlation
gws.cgwb.merged.cor = 
  gws.cgwb.merged %>%
  mutate(GWSA_mean = as.numeric(GWSA_mean),
         measured.mean = as.numeric(measured.mean)) %>%
  group_by(rowname) %>% 
  drop_na() %>%
  dplyr::summarise(cor1=cor(GWSA_mean, measured.mean))

#Link the row names to vector
gws.cgwb.cor.dug1 = 
  vectorTemplate %>%
  merge(gws.cgwb.merged.cor, by = 'rowname') %>%
  fasterize(rasterTemplate, field = 'cor1')

############################################################################################

#Load the in-situ CGWB data
cgwb.bore = 
  fread(paste0(pathIn, 'CGWB_w_specificYield.csv')) %>%
  dplyr::select(-(`May 1996`:`Nov 2002`)) %>%
  mutate(Nas = rowSums(is.na(.[,6:62]))) %>% #Keep only the rows after 2003
  dplyr::filter(Nas<(56-24) & SITE_TYPE != 'Dug Well') %>% #Filter out wells with less than 6 years of data
  st_as_sf(coords = c("LON", "LAT"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_join(vectorTemplate[,c('rowname')]) %>%
  dplyr::select(-V1, -STATE, -DISTRICT, -SITE_TYPE, -Nas)

colname.ref = colnames(cgwb.bore)

cgwb.long = 
  cgwb.bore %>%
  st_drop_geometry() %>% as.data.table() %>%
  melt(id.vars    = c("WLCODE", "rowname", "sy_final"),
       measure.vars  = colname.ref[2:58],
       variable.name = "yearmon",
       value.name    = "measured") %>%
  group_by(WLCODE) %>%
  mutate(mean_all = mean(measured, na.rm = T)) %>% 
  rowwise() %>%
  mutate(gwl.anomaly = -(measured - mean_all)) %>% #Subtract the all-time mean to obtain a anomaly; take negative to get depth below surface
  mutate(gws.anomaly = gwl.anomaly*(sy_final/100))


#Now summarize by grid cell
cgwb.grid = 
  cgwb.long %>% 
  group_by(rowname, yearmon) %>%
  dplyr::summarise(measured.mean = mean(gws.anomaly, na.rm = T)) %>%
  as.data.table() %>%
  mutate(yearmon = as.character(yearmon))


#Link grid-cell id to grace data
gws_05.SA = 
  gws_05.v %>% 
  filter(as.character(WB_NAME) == 'India') %>%
  st_join(vectorTemplate[,c('rowname')]) %>%
  st_drop_geometry() %>% as.data.table() %>%
  data.table::melt(id.vars = c("cell_id", "rowname"),
                   measure.vars = 3:212,
                   variable.name = "ymd", value.name = "GWSA_mean") %>%
  mutate(yearmon = as.yearmon(as.Date(ymd)))

#Create a merged dataset
gws.cgwb.merged = 
  gws_05.SA %>% dplyr::select(rowname, yearmon, ymd, GWSA_mean) %>%
  mutate(yearmon = as.character(yearmon),
         ymd = as.Date(ymd)) %>%
  merge(cgwb.grid, by = c('yearmon', 'rowname')) %>%
  arrange(rowname, ymd)

#Assess the correlation
gws.cgwb.merged.cor = 
  gws.cgwb.merged %>%
  mutate(GWSA_mean = as.numeric(GWSA_mean),
         measured.mean = as.numeric(measured.mean)) %>%
  group_by(rowname) %>% 
  drop_na() %>%
  dplyr::summarise(cor1=cor(GWSA_mean, measured.mean))

#Link the row names to vector
gws.cgwb.cor.bore1 = 
  vectorTemplate %>%
  merge(gws.cgwb.merged.cor, by = 'rowname') %>%
  fasterize(rasterTemplate, field = 'cor1')



############################################################################################

color_pal = c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')

plot(gws.cgwb.cor.dug1,
     breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
     col = color_pal, 
     legend.args = list(text = 'Cor (r)', side = 4, 
                        font = 2, line = 2.5, cex = 0.8)) 

ggplot(cgwb.dug) + geom_sf(size=0.5)
hist(gws.cgwb.cor.dug1)

plot(gws.cgwb.cor.bore1,
     breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
     col = color_pal, 
     legend.args = list(text = 'Cor (r)', side = 4, 
                        font = 2, line = 2.5, cex = 0.8)) 

ggplot(cgwb.bore) + geom_sf(size=0.5)
hist(gws.cgwb.cor.bore1)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

gws_mean_05 = 
  tidync(paste0(pathData, 'GWS_whymap_025deg.nc')) %>%
  hyper_tibble() %>% as.data.table()

#Create a date sequence
date.seq = 
  seq(as.Date("2003/2/1"), as.Date("2021/9/1"), "month") %>%
  as.data.table() %>%
  rownames_to_column() %>%
  .[, rowname := as.integer(rowname)]
colnames(date.seq) = c('rowname', 'date')

#Merge the dates using the rowname column
gws_mean_05 = 
  gws_mean_05 %>%
  merge(date.seq, by.x = 'time', by.y = 'rowname', all.x = T) %>%
  .[, ':='(yearmon = substr(date, 1, 7))] %>%
  dplyr::rename(GWSA_mean = GWS_whymap_025deg)

#Convert to wide
gws_05_wide = 
  gws_mean_05 %>% 
  dplyr::select(GWSA_mean, lat, lon, date) %>%
  tidyr::spread(key = date, value = GWSA_mean) %>%
  mutate(cell_id = paste0(lon, lat)) %>%
  as.data.table()

#Convert to spatial object
gws_05.v = 
  st_as_sf(gws_05_wide, coords = c("lon", "lat"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs")

#Merge country data with GRACE
gws_05.v = 
  gws_05.v %>%
  st_make_valid() %>%
  st_join(wb_regions)

#Need to isolate south asia for plots
gws_05_sa = 
  gws_05.v %>% filter(as.character(WB_NAME) == 'India') 

gws_temp = gws_05_wide %>% filter(cell_id %in% gws_05_sa$cell_id)

gws.raster = 
  rasterFromXYZ(gws_temp[,.(lon, lat, `2003-02-01`)])
#values(gws.raster) = 1


############################################################################################
#Need to essentially use the gws.raster to create a fishnet polygon
#which will be used to compare the output from GRACE and in-situ
rasterTemplate = gws.raster
vectorTemplate = 
  rasterTemplate %>%
  st_as_stars() %>%
  st_as_sf(as_points = FALSE, crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  rownames_to_column()

st_crs(vectorTemplate) = "+proj=longlat +datum=WGS84 +no_defs"

############################################################################################

#Load the in-situ CGWB data
cgwb.dug = 
  fread(paste0(pathIn, 'CGWB_w_specificYield.csv')) %>%
  dplyr::select(-(`May 1996`:`Nov 2002`)) %>%
  mutate(Nas = rowSums(is.na(.[,6:62]))) %>% #Keep only the rows after 2003
  dplyr::filter(Nas<(56-24) & SITE_TYPE == 'Dug Well') %>% #Filter out wells with less than 6 years of data
  st_as_sf(coords = c("LON", "LAT"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_join(vectorTemplate[,c('rowname')]) %>%
  dplyr::select(-V1, -STATE, -DISTRICT, -SITE_TYPE, -Nas)

colname.ref = colnames(cgwb.dug)

cgwb.long = 
  cgwb.dug %>%
  st_drop_geometry() %>% as.data.table() %>%
  melt(id.vars    = c("WLCODE", "rowname", "sy_final"),
       measure.vars  = colname.ref[2:58],
       variable.name = "yearmon",
       value.name    = "measured") %>%
  group_by(WLCODE) %>%
  mutate(mean_all = mean(measured, na.rm = T)) %>% 
  rowwise() %>%
  mutate(gwl.anomaly = -(measured - mean_all)) %>% #Subtract the all-time mean to obtain a anomaly; take negative to get depth below surface
  mutate(gws.anomaly = gwl.anomaly*(sy_final/100))


#Now summarize by grid cell
cgwb.grid = 
  cgwb.long %>% 
  group_by(rowname, yearmon) %>%
  dplyr::summarise(measured.mean = mean(gws.anomaly, na.rm = T)) %>%
  as.data.table() %>%
  mutate(yearmon = as.character(yearmon))


#Link grid-cell id to grace data
gws_05.SA = 
  gws_05.v %>% 
  filter(as.character(WB_NAME) == 'India') %>%
  st_join(vectorTemplate[,c('rowname')]) %>%
  st_drop_geometry() %>% as.data.table() %>%
  data.table::melt(id.vars = c("cell_id", "rowname"),
                   measure.vars = 3:212,
                   variable.name = "ymd", value.name = "GWSA_mean") %>%
  mutate(yearmon = as.yearmon(as.Date(ymd)))

#Create a merged dataset
gws.cgwb.merged = 
  gws_05.SA %>% dplyr::select(rowname, yearmon, ymd, GWSA_mean) %>%
  mutate(yearmon = as.character(yearmon),
         ymd = as.Date(ymd)) %>%
  merge(cgwb.grid, by = c('yearmon', 'rowname')) %>%
  arrange(rowname, ymd)

#Assess the correlation
gws.cgwb.merged.cor = 
  gws.cgwb.merged %>%
  mutate(GWSA_mean = as.numeric(GWSA_mean),
         measured.mean = as.numeric(measured.mean)) %>%
  group_by(rowname) %>% 
  drop_na() %>%
  dplyr::summarise(cor1=cor(GWSA_mean, measured.mean))

#Link the row names to vector
gws.cgwb.cor.dug2 = 
  vectorTemplate %>%
  merge(gws.cgwb.merged.cor, by = 'rowname') %>%
  fasterize(rasterTemplate, field = 'cor1')

############################################################################################

#Load the in-situ CGWB data
cgwb.bore = 
  fread(paste0(pathIn, 'CGWB_w_specificYield.csv')) %>%
  dplyr::select(-(`May 1996`:`Nov 2002`)) %>%
  mutate(Nas = rowSums(is.na(.[,6:62]))) %>% #Keep only the rows after 2003
  dplyr::filter(Nas<(56-24) & SITE_TYPE != 'Dug Well') %>% #Filter out wells with less than 6 years of data
  st_as_sf(coords = c("LON", "LAT"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_join(vectorTemplate[,c('rowname')]) %>%
  dplyr::select(-V1, -STATE, -DISTRICT, -SITE_TYPE, -Nas)

colname.ref = colnames(cgwb.bore)

cgwb.long = 
  cgwb.bore %>%
  st_drop_geometry() %>% as.data.table() %>%
  melt(id.vars    = c("WLCODE", "rowname", "sy_final"),
       measure.vars  = colname.ref[2:58],
       variable.name = "yearmon",
       value.name    = "measured") %>%
  group_by(WLCODE) %>%
  mutate(mean_all = mean(measured, na.rm = T)) %>% 
  rowwise() %>%
  mutate(gwl.anomaly = -(measured - mean_all)) %>% #Subtract the all-time mean to obtain a anomaly; take negative to get depth below surface
  mutate(gws.anomaly = gwl.anomaly*(sy_final/100))


#Now summarize by grid cell
cgwb.grid = 
  cgwb.long %>% 
  group_by(rowname, yearmon) %>%
  dplyr::summarise(measured.mean = mean(gws.anomaly, na.rm = T)) %>%
  as.data.table() %>%
  mutate(yearmon = as.character(yearmon))


#Link grid-cell id to grace data
gws_05.SA = 
  gws_05.v %>% 
  filter(as.character(WB_NAME) == 'India') %>%
  st_join(vectorTemplate[,c('rowname')]) %>%
  st_drop_geometry() %>% as.data.table() %>%
  data.table::melt(id.vars = c("cell_id", "rowname"),
                   measure.vars = 3:212,
                   variable.name = "ymd", value.name = "GWSA_mean") %>%
  mutate(yearmon = as.yearmon(as.Date(ymd)))

#Create a merged dataset
gws.cgwb.merged = 
  gws_05.SA %>% dplyr::select(rowname, yearmon, ymd, GWSA_mean) %>%
  mutate(yearmon = as.character(yearmon),
         ymd = as.Date(ymd)) %>%
  merge(cgwb.grid, by = c('yearmon', 'rowname')) %>%
  arrange(rowname, ymd)

#Assess the correlation
gws.cgwb.merged.cor = 
  gws.cgwb.merged %>%
  mutate(GWSA_mean = as.numeric(GWSA_mean),
         measured.mean = as.numeric(measured.mean)) %>%
  group_by(rowname) %>% 
  drop_na() %>%
  dplyr::summarise(cor1=cor(GWSA_mean, measured.mean))

#Link the row names to vector
gws.cgwb.cor.bore2 = 
  vectorTemplate %>%
  merge(gws.cgwb.merged.cor, by = 'rowname') %>%
  fasterize(rasterTemplate, field = 'cor1')




############################################################################################

color_pal = c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')

plot(gws.cgwb.cor.dug2,
     breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
     col = color_pal, 
     legend.args = list(text = 'Cor (r)', side = 4, 
                        font = 2, line = 2.5, cex = 0.8)) 

ggplot(cgwb.dug) + geom_sf(size=0.5)
hist(gws.cgwb.cor.dug2)
summary(gws.cgwb.cor.dug2)

plot(gws.cgwb.cor.bore2,
     breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
     col = color_pal, 
     legend.args = list(text = 'Cor (r)', side = 4, 
                        font = 2, line = 2.5, cex = 0.8)) 

ggplot(cgwb.bore) + geom_sf(size=0.5)
hist(gws.cgwb.cor.bore2)
summary(gws.cgwb.cor.bore2)

############################################################################################
writeRaster(gws.cgwb.cor.dug1, paste0(pathOut, 'Cor_WB_dug.tif'), overwrite=TRUE)
writeRaster(gws.cgwb.cor.bore1, paste0(pathOut, 'Cor_WB_bore.tif'), overwrite=TRUE)
writeRaster(gws.cgwb.cor.dug2, paste0(pathOut, 'Cor_whymap_dug.tif'), overwrite=TRUE)
writeRaster(gws.cgwb.cor.bore2, paste0(pathOut, 'Cor_whymap_bore.tif'), overwrite=TRUE)

