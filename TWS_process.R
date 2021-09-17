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

tws_lwe_base = cbind(tws_lwe[,1:3], tws_lwe_base)
  
#fwrite(tws_lwe_base, "GRACE_TWS_scaled.csv")

tws_trends = 
  apply(tws_lwe_base, 1, zyp.yuepilon)

tws_trends_estimates = 
  tws_lwe[,1:3] %>%
  mutate(trends = tws_trends[2,], sig = tws_trends[6,]) %>%
  mutate(lon2 = ifelse(lon>180, lon-360, lon))

#fwrite(tws_trends_estimates, "GRACE_TWS_trends.csv")


test = 
  tidync('/Users/tejasvi/Dropbox/Mac/Downloads/GLDAS/GLDAS_NOAH025_M.A200309.021.nc4.SUB.nc4')


gldas_noah = 
  tidync('~/Dropbox/Mac/Downloads/GLDAS/gldas_noah_merge.nc4') %>%
  hyper_filter(lon = lon > 65 & lon < 100,  lat = lat > 5 & lat < 40, time = time < 360) %>%
  hyper_tibble() %>% as.data.table() 

#fwrite(test, "GLDAS.csv")



