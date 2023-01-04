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

#Open the reference
ref.tws = fread(paste0(proj_dir, "GRACE_Data/Reference/GRACE_JPL_TWS_1degree_BSL2017_220725.csv"))

#The GRACE Mascon dataset
tws_mascon =
  tidync(paste0(proj_dir,'GRACE_data/GSFC_Mascons/GSFC_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02_1degree.nc')) %>%
  activate("D2,D3,D0") %>%
  hyper_tibble() %>% as.data.table() %>%
  .[order(lat, lon)] %>% #ensures the order is the same
  group_by(lat,lon) %>%
  mutate(cell_id = cur_group_id())  %>% as.data.table()

# dates don't align with the JPL solutions by 2 days
tws_mascon[time == 6148,]$time = 6148 - 2

#Dates 
dates1 = 
  tws_mascon[,.(time)]%>% distinct() %>%
  mutate(dates = as.Date(time, origin = '2002-01-01')) %>%
  mutate(ym = substr(ymd(dates), 1, 7)) 

dates1 = dates1[!duplicated(dates1$ym)]

########################
######### Analysis
########################
#Add dates to the TWS df; also scale the values 
tws_mascon = 
  tws_mascon %>% 
  merge(dates1[,.(ym, time)], by = 'time') %>%
  .[, lwe_corr := lwe_thickness] 

#Convert from long to wide format
tws_lwe = 
  tws_mascon[,.(lwe_corr, lon, lat, cell_id, ym)] %>%
  dcast(lon + lat + cell_id ~ ym, value.var = "lwe_corr") 

#Check if colnames match reference JPL dataset
colnames(ref.tws) %in% colnames(tws_lwe)

#Estimate the new baseline (2002-2020)
bsl_colnum = 166 #199 = 2020; 165 = 2017; 118 = 2012
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
       paste0(proj_dir, "GRACE_Data/GSFC_Mascons/GRACE_GSFC_TWS_1degree_221230.csv"))

#Now add this dataset to the reference sheet (FileSummary.csv)
