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

#Can pick the z-score or cm-equivalent version of GWS

z_score = F

if (z_score == T){
  #Load the GWS data
  gws_TS = 
    fread(paste0(proj_dir, 
                 'Outputs/GWS/9003_GRACE_GWS_TWS_3Mascons_SWS_Ensemble_2002_2020_BSL2017.csv'))
  
  #Estimate the z-scores
  gws_TS_mean =  apply(gws_TS[, 4:193], 1, mean) 
  gws_TS_sd = apply(gws_TS[, 4:193], 1, sd) 
  
  #Add the columns for lat/long/cell_id
  gws_TS_zscore = 
    (gws_TS[, 4:193] - gws_TS_mean)/gws_TS_sd
  
  gws_TS_zscore = 
    cbind(gws_TS[,1:3], gws_TS_zscore) 
  
  grace = gws_TS_zscore
  
} else {
  grace = 
    fread(paste0(proj_dir, 
                 'Outputs/GWS/9003_GRACE_GWS_TWS_3Mascons_SWS_Ensemble_2002_2020_BSL2017.csv'))
}

#Source: https://doi.org/10.1038/s43016-021-00429-z
crop_area = 
  raster("/Users/tejasvi/Dropbox/WB/Cropped Area/Global_cropland_3km_2019.tif")

grace_grid = 
  grace %>% 
  dplyr::select(1:3) %>%
  rasterFromXYZ(crs = crs(crop_area)) 

#Aggregate the cropped area grid to the grace grid using billinear interpolation
#Convert to pixel data frame so it can be merged with GRACE grid
crop_area_grace = 
  crop_area %>%
  resample(grace_grid) %>%
  rasterToPoints() %>%
  as.data.table() %>%
  dplyr::rename(Per_crop_area = Global_cropland_3km_2019)

grace_wcrop = merge(grace, crop_area_grace, 
                    by.x = c('lon', 'lat'), by.y = c('x','y'), all.x = T)

#Make sure the merging was proper -- regions with expected high cropped area 
# show up (e.g. India)
grace_test = 
  grace_wcrop %>% dplyr::select(lon, lat, Per_crop_area) %>%
  rasterFromXYZ(crs = crs(crop_area)) %>%
  plot()

####################
#Use the country level median estimating script
# 
# #Load the data --this time it contains %cropped area (z-score version)
# gws_TS = 
#   fread(paste0(pathOut,'GRACE_GWS_2002_2020_zscore.csv'))

gws_TS = grace_wcrop
######################################################

#Convert to spatial object
gws_TS_spatial = 
  st_as_sf(gws_TS, coords = c("lon", "lat"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs")

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

crop_flag = T
if(crop_flag == T) {
  gws_TS_country =
    gws_TS_country %>%
    dplyr::filter(Per_crop_area>20)
}

#Test to see how the remaining GRACE points look
plot(gws_TS_country['Per_crop_area'], cex = 0.25)

#Get list of countries with atleast 36 points
ag.countries = 
  gws_TS_country %>%
  as.data.frame() %>%
  group_by(WB_NAME) %>%
  summarise(count = n()) %>%
  filter(count > 36)

gws.country.50 = 
  gws_TS_country %>%
  dplyr::mutate(area_flag = 
                  ifelse(WB_NAME %in% ag.countries$WB_NAME, T, F)) %>%
  dplyr::select(-Per_crop_area) %>%
  group_by(WB_NAME) %>%
  summarise_if(is.numeric, median, na.rm = TRUE)

gws.country.25 = 
  gws_TS_country %>%
  dplyr::mutate(area_flag = 
                  ifelse(WB_NAME %in% ag.countries$WB_NAME, T, F)) %>%
  dplyr::select(-Per_crop_area) %>%
  group_by(WB_NAME) %>%
  summarise_if(is.numeric, function (x){quantile(x,probs = 0.25, na.rm = TRUE)})

gws.country.75 = 
  gws_TS_country %>%
  dplyr::mutate(area_flag = 
                  ifelse(WB_NAME %in% ag.countries$WB_NAME, T, F)) %>%
  dplyr::select(-Per_crop_area) %>%
  group_by(WB_NAME) %>%
  summarise_if(is.numeric, function (x){quantile(x,probs = 0.75, na.rm = TRUE)})

gws.50.long = 
  gws.country.50 %>%
  gather(yearmon, gws_median, `2002-04`:`2020-12`) %>%
  as.data.table() %>%
  dplyr::select(-geometry, -`2021-01`) %>%
  drop_na()


gws.25.long = 
  gws.country.25 %>%
  gather(yearmon, gws_25, `2002-04`:`2020-12`) %>%
  as.data.table() %>%
  dplyr::select(-geometry, -`2021-01`) %>%
  drop_na()


gws.75.long = 
  gws.country.75 %>%
  gather(yearmon, gws_75, `2002-04`:`2020-12`) %>%
  as.data.table() %>%
  dplyr::select(-geometry, -`2021-01`) %>%
  drop_na()

gws.comb = 
  merge(gws.50.long, gws.25.long, by = c('WB_NAME','yearmon'), all.x = T) %>%
  merge(gws.75.long, by = c('WB_NAME','yearmon'), all.x = T) %>%
  merge(wb_regions %>% as.data.frame() %>% dplyr::select(-geometry), 
        by = 'WB_NAME', all.x = T) %>%
  dplyr::select(-TYPE) %>%
  dplyr::mutate(area_flag = 
                  ifelse(WB_NAME %in% ag.countries$WB_NAME, T, F)) 


gws.wide.country = 
  gws.comb %>%
  dplyr::select(1:3) %>%
  spread(WB_NAME, gws_median) 

gws.xts = 
  gws.wide.country[, 2:ncol(gws.wide.country)] %>%
  #  dplyr::select(-yearmon) %>%
  xts(as.yearmon(gws.wide.country$yearmon))

gws.smooth = 
  rollmean(gws.xts, k = 24) %>%
  as.data.table()


# fwrite(gws.wide.country, 'Country_level_TS/country_level_gws_MEDIAN.csv')
# fwrite(gws.smooth, 'Country_level_TS/country_level_gws_MEDIAN_smooth.csv')
# fwrite(gws.comb, 'Output/cm_equivalent/country_level_gws_centimeter_crop.csv')

#Name Creation
if(z_score == T){
  name1 = 'zscore'
} else{name1 = 'centimeter'}

if(crop_flag == T){
  name2 = 'crop'
} else{name2 = ''}

outName = paste('country_level_gws', name1, name2, sep = '_')

fwrite(gws.comb, 
       paste0(pathOut, outName, '.csv'))



plot(gws.xts$Canada, add = T)
plot(test$Canada)






test = 
  fread('Country_level_TS/country_level_gws_MEDIAN.csv') %>%
  dplyr::select(2:162) %>%
  xts(as.yearmon(gws.wide.country$yearmon))


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

