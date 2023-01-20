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
pathOut =  "/Users/tejasvi/Dropbox/WB/GRACE-Subsidies/"
pathIn = '/Users/tejasvi/Dropbox/gwflagship_typologies/'

############################################################################################
#Load the GWS data
gws_TS = 
  fread(paste0(proj_dir, 
               'Outputs/GWS/9003_GRACE_GWS_TWS_3Mascons_SWS_Ensemble_2002_2020_BSL2017.csv'))

#Estimate the z-scores
gws_TS_mean =  apply(gws_TS[, 4:164], 1, mean) 
gws_TS_sd = apply(gws_TS[, 4:164], 1, sd) 

#Add the columns for lat/long/cell_id
gws_TS_zscore = 
  (gws_TS[, 4:193] - gws_TS_mean)/gws_TS_sd

gws_TS_zscore = 
  cbind(gws_TS[,1:3], gws_TS_zscore) 

#fwrite(gws_TS_zscore, paste0(pathOut,'GRACE_GWS_2002_2017_zscore.csv'))


#Convert to spatial object
gws_TS_spatial = 
  st_as_sf(gws_TS_zscore, coords = c("lon", "lat"), 
           crs = 4326)

####Load World Regions
wb_regions = 
  st_read(paste0(proj_dir, "Spatial Files/WB_Regions/WB_countries_Admin0_10m.shp")) %>%
  dplyr::select(WB_NAME, ISO_A2, ISO_A3, ISO_N3, TYPE, REGION_WB) %>%
  filter(TYPE != 'Dependency') %>%
  st_make_valid()

wb_regions_ns = 
  wb_regions %>% as.data.table() %>% dplyr::select(-geometry) %>% distinct()

#Merge country data with GRACE
gws_TS_country = 
  gws_TS_spatial %>%
  st_make_valid() %>%
  st_join(wb_regions)

gws.country.50 = 
  gws_TS_country %>%
  group_by(WB_NAME) %>%
  summarise_if(is.numeric, median, na.rm = TRUE)

gws.country.25 = 
  gws_TS_country %>%
  group_by(WB_NAME) %>%
  summarise_if(is.numeric, function (x){quantile(x,probs = 0.25, na.rm = TRUE)})

gws.country.75 = 
  gws_TS_country %>%
  group_by(WB_NAME) %>%
  summarise_if(is.numeric, function (x){quantile(x,probs = 0.75, na.rm = TRUE)})

gws.50.long = 
  gws.country.50 %>%
  gather(yearmon, gws_median, `2002-04`:`2020-12`) %>%
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  drop_na()
  
gws.25.long = 
  gws.country.25 %>%
  gather(yearmon, gws_25, `2002-04`:`2020-12`) %>%
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  drop_na()

gws.75.long = 
  gws.country.75 %>%
  gather(yearmon, gws_75, `2002-04`:`2020-12`) %>%
  as.data.table() %>%
  dplyr::select(-geometry) %>%
  drop_na()

gws.comb = 
  merge(gws.50.long, gws.25.long, by = c('WB_NAME','yearmon'), all.x = T) %>%
  merge(gws.75.long, by = c('WB_NAME','yearmon'), all.x = T) %>%
  merge(wb_regions %>% as.data.frame() %>% dplyr::select(-geometry), 
        by = 'WB_NAME', all.x = T) %>%
  dplyr::select(-TYPE)

gws.wide.country = 
  gws.comb %>%
  dplyr::select(1:3) %>%
  spread(WB_NAME, gws_median) 

gws.xts = 
  gws.wide.country[, 2:154] %>%
#  dplyr::select(-yearmon) %>%
  xts(as.yearmon(gws.wide.country$yearmon))

gws.smooth = 
  rollmean(gws.xts, k = 24) %>%
  as.data.table()
  

# fwrite(gws.wide.country, 'Country_level_TS/country_level_gws_MEDIAN.csv')
# fwrite(gws.smooth, 'Country_level_TS/country_level_gws_MEDIAN_smooth.csv')
fwrite(gws.comb, paste0(pathOut, 'Country_level_GWS_zscore_230116.csv'))


#India test

india.st = 
  st_read("States/Admin2.shp") %>%
  st_make_valid()

test = 
  gws_TS_country %>% 
  filter(WB_NAME=='India') %>%
  st_join(india.st) %>%
  filter(ST_NM=='Punjab')

test.50 = 
  test %>%
  summarise_if(is.numeric, median, na.rm = TRUE) %>%
  as.data.frame()

plot(t(test.50[,3:160]))
