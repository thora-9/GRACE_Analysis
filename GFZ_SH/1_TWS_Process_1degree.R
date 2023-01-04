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

proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"

#Load the datasets

#The GRACE SH dataset
tws_mascon1 =
  tidync(paste0(proj_dir,'GRACE_data/GFZ_Solutions/GeodeticsGravity_TELLUS_GRAC_L3_GFZ_RL06_LND_v04_CORR_GRID.nc')) %>%
  activate("D1,D2,D0") %>%
  hyper_tibble() %>% as.data.table() %>%
  .[order(lat, lon)] %>% #ensures the order is the same
  group_by(lat,lon) %>%
  mutate(cell_id = cur_group_id())  %>% as.data.table()


#Dates 
dates1a = 
  tws_mascon1[,.(time)] %>% distinct() %>%
  mutate(time2 = time)

#2 of the dates don't align with the JPL solutions by 1 day $ 4 days
#Manually add those days to make them align

dates1a[111,]$time = dates1a[111,]$time + 1
dates1a[112,]$time = dates1a[112,]$time + 4


dates1a = 
  dates1a %>%
  mutate(dates = as.Date(time, origin = '2002-01-01')) %>%
  mutate(ym = substr(ymd(dates), 1, 7)) 

dates1a = dates1a[!duplicated(dates1a$ym)]

#The GRACE SH dataset
tws_mascon2 =
  tidync(paste0(proj_dir,'GRACE_data/GFZ_Solutions/GeodeticsGravity_TELLUS_GRFO_L3_GFZ_RL06_LND_v04_CORR_GRID.nc')) %>%
  activate("D1,D2,D0") %>%
  hyper_tibble() %>% as.data.table() %>%
  .[order(lat, lon)] %>% #ensures the order is the same
  group_by(lat,lon) %>%
  mutate(cell_id = cur_group_id())  %>% as.data.table()


#Dates 
dates1b = 
  tws_mascon2[,.(time)] %>% distinct() %>% 
  mutate(time2 = time) %>%
  mutate(dates = as.Date(time, origin = '2002-01-01')) %>%
  mutate(ym = substr(ymd(dates), 1, 7)) 

dates1b = dates1b[!duplicated(dates1b$ym)]


#Merge the two
tws_mascon = 
  rbind(tws_mascon1, tws_mascon2) %>%
  dplyr::select(-cell_id) %>%
  group_by(lat,lon) %>%
  mutate(cell_id = cur_group_id())  %>% as.data.table() 

dates1 = rbind(dates1a, dates1b)

########################
######### Analysis
########################
#Add dates to the TWS df; also scale the values 
#Convert units from 'm' to 'cm'
tws_mascon = 
  tws_mascon %>% merge(dates1[,.(ym, time = time2)], by = 'time') %>%
  .[, lwe_corr := lwe_thickness*100] 

#Convert from long to wide format
tws_lwe = 
  tws_mascon[,.(lwe_corr, lon, lat, cell_id, ym)] %>%
  dcast(lon + lat + cell_id ~ ym, value.var = "lwe_corr") 

#Estimate the new baseline (2002-2020)
bsl_colnum = 164 #199 = 2020; 165 = 2012; 118 = 2012
mean_tws = apply(tws_lwe[, 4:bsl_colnum], 1, mean) 

#Estimate the new baseline
sd_tws = apply(tws_lwe[, 4:bsl_colnum], 1, sd) 

#Remove the new baseline (currently not using data from GRACE-FO)
tws_lwe_base = 
  tws_lwe[, 4:199] - mean_tws 

#Estimate the new baseline
#sd_tws2 = apply(tws_lwe_base, 1, sd) 

tws_lwe_base = 
  cbind(tws_lwe[,1:3], tws_lwe_base) %>%
  mutate(lon = ifelse(lon>180, lon-360, lon))

fwrite(tws_lwe_base,
       paste0(proj_dir, "GRACE_data/GFZ_Solutions/GRACE_GFZ_TWS_1degree_220821.csv"))

