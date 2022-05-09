library(tidync)
library(data.table)
library(tidyverse)
library(lubridate)
library(zyp)

tws_mascon =
  tidync('~/Dropbox/WB/GRACE_Analysis/GRCTellus.JPL.200204_202106.GLO.RL06M.MSCNv02CRI.nc') %>%
  activate("D0,D1,D2") %>%
  hyper_tibble() %>% as.data.table() %>%
  .[order(lat, lon)] %>% #ensures the order is the same
  group_by(lat,lon) %>%
  mutate(cell_id = cur_group_id())  %>% as.data.table()

mascons = 
  tidync('~/Dropbox/WB/GRACE_Analysis/JPL_MSCNv02_PLACEMENT.nc') %>%
  activate("D1") %>%
  hyper_tibble() %>% as.data.table()

clm_factors = 
  tidync('~/Dropbox/WB/GRACE_Analysis/CLM4.SCALE_FACTOR.JPL.MSCNv02CRI.nc') %>%
  activate("D0,D1") %>%
  hyper_tibble() %>% as.data.table()

land_filter = 
  tidync('~/Dropbox/WB/GRACE_Analysis/LAND_MASK.CRI.nc') %>%
  activate("D0,D1") %>%
  hyper_tibble() %>% as.data.table() 

dates1 = 
  tws_mascon[,.(time)]%>% distinct() %>%
  mutate(dates = as.Date(time, origin = '2002-01-01')) %>%
  mutate(ym = substr(ymd(dates), 1, 7)) 

dates1 = dates1[!duplicated(dates1$ym)]

tws_mascon = 
  tws_mascon %>% merge(dates1[,.(ym, time)], by = 'time') %>%
  merge(clm_factors, by = c("lat", "lon")) %>% merge(land_filter, by = c("lat", "lon")) %>%
  .[, lwe_corr := scale_factor*lwe_thickness] 
  
tws_lwe = 
  tws_mascon[,.(lwe_corr, lon, lat, cell_id, ym)] %>%
  dcast(lon + lat + cell_id ~ ym, value.var = "lwe_corr") 

#Estimate the new baseline
mean_tws = apply(tws_lwe[, 4:164], 1, mean) 

#Remove the new baseline (currently not using data from GRACE-FO)
tws_lwe_base = 
  tws_lwe[, 4:164] - mean_tws 

tws_lwe_base = 
  cbind(tws_lwe[,1:3], tws_lwe_base) %>%
  mutate(lon = ifelse(lon>180, lon-360, lon))
  
#fwrite(tws_lwe_base, "GRACE_TWS_scaled.csv")

tws_trends = 
  apply(tws_lwe_base[, 3:ncol(tws_lwe_base)], 1, zyp.yuepilon)

tws_trends_estimates = 
  tws_lwe[,1:3] %>%
  mutate(trends = tws_trends[2,], sig = tws_trends[6,]) %>%
  mutate(lon2 = ifelse(lon>180, lon-360, lon))

#fwrite(tws_trends_estimates, "GRACE_TWS_trends.csv")


test = 
  tidync('/Users/tejasvi/Dropbox/Mac/Downloads/GLDAS/GLDAS_NOAH025_M.A200309.021.nc4.SUB.nc4')


#Use the methodology developed by Shamsuddhua and Taylor (2020)

# For accumulated variables such as Qs_acc, the monthly mean surface runoff is the 
# average 3-hour accumulation over all 3-hour intervals in April 1979. To compute 
# monthly accumulation, use this formula: 
#   
#   Qs_acc (April){kg/m2} = Qs_acc (April){kg/m2/3hr} * 8{3hr/day} * 30{days} 

gldas_noah = 
  tidync('~/Dropbox/Mac/Downloads/GLDAS/gldas_noah_merge.nc4') %>%
  hyper_tibble() %>% as.data.table() %>%
  .[,Qs_month := Qs_acc * 8 * 30] %>% #Convert accumulated runoff into monthly runoff (see formula above)
  .[, total_water_without_runoff := 
      (SWE_inst + SoilMoi0_10cm_inst + SoilMoi10_40cm_inst +
      SoilMoi40_100cm_inst + SoilMoi100_200cm_inst + CanopInt_inst)/10] %>%
  .[, total_water := 
      (Qs_month + SWE_inst + SoilMoi0_10cm_inst + SoilMoi10_40cm_inst +
      SoilMoi40_100cm_inst + SoilMoi100_200cm_inst + CanopInt_inst)/10] 
# convert mass to volume: vol = kg / (1000 kg/m3) [Density of water] -> vol = (1/1000) m3
# Therefore: kg/m2 -> (1/1000) [m3/m2] -> (1/10) [cm]

gldas_dates =
  copy(gldas_noah[,.(time)]) %>% distinct() %>%
  .[,ym := substr(as.Date(time, origin = '2000-01-01'), 1, 7)]  

gldas_noah = 
  gldas_noah %>% merge(gldas_dates, by = 'time') %>%
  .[ym %in% dates1$ym] 

##########Without runoff
gldas_noah_spread = 
  gldas_noah[,.(total_water_without_runoff, lon, lat, ym)] %>%
  dcast(lon + lat ~ ym, value.var = "total_water_without_runoff") 

#Estimate the new baseline land water content
mean_gldas = apply(gldas_noah_spread[, 3:163], 1, mean) 

#Remove the new baseline (currently not using data from GRACE-FO)
gldas_lwc_base = 
  gldas_noah_spread[, 3:163] - mean_gldas

gldas_noah_base_rem = 
  cbind(gldas_noah_spread[,1:2], gldas_lwc_base) %>%
  merge(tws_lwe_base[,.(lon,lat,cell_id)], by = c("lat", "lon"), all.x = F)

#fwrite(gldas_noah_base_rem, "GLDAS_2002_2017_SWS_without_runoff.csv")

###########With runoff
gldas_noah_spread = 
  gldas_noah[,.(total_water, lon, lat, ym)] %>%
  dcast(lon + lat ~ ym, value.var = "total_water") 

#Estimate the new baseline land water content
mean_gldas = apply(gldas_noah_spread[, 3:163], 1, mean) 

#Remove the new baseline (currently not using data from GRACE-FO)
gldas_lwc_base = 
  gldas_noah_spread[, 3:163] - mean_gldas

gldas_noah_base_rem = 
  cbind(gldas_noah_spread[,1:2], gldas_lwc_base) %>%
  merge(tws_lwe_base[,.(lon,lat,cell_id)], by = c("lat", "lon"), all.x = F)

#fwrite(gldas_noah_base_rem, "GLDAS_2002_2017_SWS_with_runoff.csv")

###Create Groundwater Storage Anomaly

sws_anomaly = 
  gldas_noah_base_rem %>% .[order(lon,lat)]

tws_anomaly = 
  tws_lwe_base[cell_id %in% sws_anomaly$cell_id] %>% .[order(lon,lat)]

identical(sws_anomaly$cell_id, tws_anomaly$cell_id)
identical(colnames(sws_anomaly)[3:163], colnames(tws_anomaly)[4:164])
identical(dim(tws_anomaly[ ,4:164]), dim(sws_anomaly[ ,3:163]))

gws_anomaly = 
  tws_anomaly[ ,4:164] - sws_anomaly[ ,3:163]

gws_anomaly = cbind(tws_anomaly[,1:3], gws_anomaly)

#fwrite(gws_anomaly, "GWS_2002_2017_with_runoff.csv")

####Trend Analysis
gws_trends = 
  apply(gws_anomaly[, 4:ncol(tws_lwe_base)], 1, zyp.yuepilon)

gws_trends_estimates = 
  gws_anomaly[,1:3] %>%
  mutate(trends = gws_trends[2,], sig = gws_trends[6,]) %>%
  mutate(lon = ifelse(lon>180, lon-360, lon))

#fwrite(gws_trends_estimates, "GRACE_GWS_trends.csv")


