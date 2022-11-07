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

proj_dir = "/Users/tejasvi/Dropbox/WB/GRACE_Ensemble/"

pathData = '/Users/tejasvi/Dropbox/gwflagship_GRACEdownscaling/Downscaled TWS_GWS v2/'

####Load World Regions
wb_regions = 
  st_read(paste0(proj_dir, "Spatial Files/WB_Regions/WB_countries_Admin0_10m.shp")) %>%
  dplyr::select(WB_NAME, ISO_A2, ISO_A3, ISO_N3, TYPE, REGION_WB) %>%
  filter(TYPE != 'Dependency') %>%
  st_make_valid()

#load the non-downscaled GRACE data
gws.OG = 
  fread(paste0(proj_dir, 
               'Outputs/GWS/GRACE_GWS_TWS_Ensemble_SWS_Ensemble_2002_2020_BSL2017.csv')) %>%
  #mutate(lon = lon + 0.001, lat = lat + 0.001) %>% #Jitter so points fall within a single polygon
  st_as_sf(coords = c("lon", "lat"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_make_valid() %>%
  st_join(wb_regions) %>%
  #st_join(vectorTemplate[,c('rowname')]) %>%
  drop_na()

#Land Filter
land.filter = 
  tidync('~/Dropbox/WB/GRACE_Analysis/LAND_MASK.CRI.nc') %>%
  activate("D0,D1") %>%
  hyper_tibble() %>% as.data.table() %>%
  mutate(lon = ifelse(lon>180, lon-360, lon))

land.filter.r = 
  rasterFromXYZ(land.filter[,.(lon, lat, land_mask)])
crs(land.filter.r) = crs(wb_regions)

#Set non-land values to NA
land.filter.r[values(land.filter.r)==0] = NA

#Convert to stars and then sf object to allow masking
land.filter.v = 
  st_as_stars(land.filter.r) %>%
  st_as_sf()

#land.filter.v.dissolved = st_union(land.filter.v)

#3 degree grid
tws.3 = 
  tidync(paste0('/Users/tejasvi/Dropbox/gwflagship_GRACEdownscaling/Downscaled TWS/', 
                'TWS_3deg.nc')) %>%
  hyper_tibble() %>% as.data.table() %>%
  dplyr::select(lat, lon) %>%
  unique() 

tws.3.r = 
  rasterFromXYZ(tws.3[,.(lon, lat)])
values(tws.3.r) = runif(5640) 
crs(tws.3.r) = crs(wb_regions)

#Crop the land.filter raster to match extents 
land.filter.3 = projectRaster(land.filter.r, tws.3.r)

pathOut =  "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"

writeRaster(land.filter.3, paste0(pathOut, '3degree_placement.tif'),
            format="GTiff", overwrite=TRUE)

