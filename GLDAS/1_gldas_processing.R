##############################
##Process GLDAS Data
##############################

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

proj_dir = "~/Dropbox/WB/GRACE_Ensemble/GLDAS_data/"

#Reference needs to be the same resolution as the GLDAS
#Load the GRACE dataset as a base
grace_ref = 
  fread(paste0("~/Dropbox/WB/GRACE_Ensemble/GRACE_Data/Reference/GRACE_JPL_TWS_1degree_BSL2017_220725.csv"))

grace_dates = 
  data.table(ym = colnames(grace_ref)[4:199])


#Use the methodology developed by Shamsuddhua and Taylor (2020)

# For accumulated variables such as Qs_acc, the monthly mean surface runoff is the 
# average 3-hour accumulation over all 3-hour intervals in April 1979. To compute 
# monthly accumulation, use this formula: 
#   
#   Qs_acc (April){kg/m2} = Qs_acc (April){kg/m2/3hr} * 8{3hr/day} * 30{days} 


gldas_noah = 
  tidync(paste0(proj_dir, 'NOAH/GLDAS_NOAH10_2000_2021.nc4')) %>%
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
  .[ym %in% grace_dates$ym] 


###########With runoff
gldas_noah_spread = 
  gldas_noah[,.(total_water, lon, lat, ym)] %>%
  dcast(lon + lat ~ ym, value.var = "total_water") 

#Estimate the new baseline (2002-2020)
bsl_colnum = 164 #199 = 2020; 165 = 2017; 118 = 2012

#Estimate the new baseline land water content
mean_gldas = apply(gldas_noah_spread[, 3:(bsl_colnum-1)], 1, mean) 

#Remove the new baseline (currently not using data from GRACE-FO)
gldas_lwc_base = 
  gldas_noah_spread[, 3:198] - mean_gldas

gldas_noah_base_rem = 
  cbind(gldas_noah_spread[,1:2], gldas_lwc_base) %>%
  merge(grace_ref[,.(lon,lat,cell_id)], by = c("lat", "lon"), all.x = F)

#Test the grid by plotting it
test_grid = 
  gldas_noah_base_rem[,c(2,1,3)] %>% 
  rasterFromXYZ(crs = crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) 

# fwrite(gldas_noah_base_rem,
#        paste0(proj_dir, 'NOAH/GLDAS_NOAH_1degree_BSL2017_SWS_with_runoff_220508.csv'))


#########################################################################################################
#########################################################################################################
#########################################################################################################
#VIC
#Use the methodology developed by Shamsuddhua and Taylor (2020)

# For accumulated variables such as Qs_acc, the monthly mean surface runoff is the 
# average 3-hour accumulation over all 3-hour intervals in April 1979. To compute 
# monthly accumulation, use this formula: 
#   
#   Qs_acc (April){kg/m2} = Qs_acc (April){kg/m2/3hr} * 8{3hr/day} * 30{days} 


gldas_vic = 
  tidync(paste0(proj_dir, 'VIC/GLDAS_VIC10_2000_2021.nc4')) %>%
  hyper_tibble() %>% as.data.table() %>%
  .[,Qs_month := Qs_acc * 8 * 30] %>% #Convert accumulated runoff into monthly runoff (see formula above)
  .[, total_water_without_runoff := 
      (SWE_inst + SoilMoi0_30cm_inst + SoilMoi_depth2_inst +
         SoilMoi_depth3_inst + CanopInt_inst)/10] %>%
  .[, total_water := 
      (Qs_month + SWE_inst + SoilMoi0_30cm_inst + SoilMoi_depth2_inst +
         SoilMoi_depth3_inst + CanopInt_inst)/10] %>%
  .[, soilM_sum := SoilMoi0_30cm_inst + SoilMoi_depth2_inst +
      SoilMoi_depth3_inst]
# convert mass to volume: vol = kg / (1000 kg/m3) [Density of water] -> vol = (1/1000) m3
# Therefore: kg/m2 -> (1/1000) [m3/m2] -> (1/10) [cm]


gldas_dates =
  copy(gldas_vic[,.(time)]) %>% distinct() %>%
  .[,ym := substr(as.Date(time, origin = '2000-01-01'), 1, 7)]  

gldas_vic = 
  gldas_vic %>% merge(gldas_dates, by = 'time') %>%
  .[ym %in% grace_dates$ym] 


###########With runoff
gldas_vic_spread = 
  gldas_vic[,.(total_water, lon, lat, ym)] %>%
  dcast(lon + lat ~ ym, value.var = "total_water") 

#Estimate the new baseline (2002-2020)
bsl_colnum = 164 #199 = 2020; 165 = 2017; 118 = 2012

#Estimate the new baseline land water content
mean_gldas = apply(gldas_vic_spread[, 3:(bsl_colnum-1)], 1, mean) 

#Remove the new baseline (currently not using data from GRACE-FO)
gldas_lwc_base = 
  gldas_vic_spread[, 3:198] - mean_gldas

gldas_vic_base_rem = 
  cbind(gldas_vic_spread[,1:2], gldas_lwc_base) %>%
  merge(grace_ref[,.(lon,lat,cell_id)], by = c("lat", "lon"), all.x = F)

#Test the grid by plotting it
test_grid = 
  gldas_vic_base_rem[,c(2,1,3)] %>% 
  rasterFromXYZ(crs = crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) 
plot(test_grid)

fwrite(gldas_vic_base_rem,
       paste0(proj_dir, 'VIC/GLDAS_VIC_1degree_BSL2017_SWS_with_runoff_220725.csv'))

#########################################################################################################
#########################################################################################################
#########################################################################################################
#CLSM
#Use the methodology developed by Shamsuddhua and Taylor (2020)

# For accumulated variables such as Qs_acc, the monthly mean surface runoff is the 
# average 3-hour accumulation over all 3-hour intervals in April 1979. To compute 
# monthly accumulation, use this formula: 
#   
#   Qs_acc (April){kg/m2} = Qs_acc (April){kg/m2/3hr} * 8{3hr/day} * 30{days} 


gldas_clsm = 
  tidync(paste0(proj_dir, 'CLSM/GLDAS_CLSM10_2000_2021.nc4')) %>%
  hyper_tibble() %>% as.data.table() %>%
  .[,Qs_month := Qs_acc * 8 * 30] %>% #Convert accumulated runoff into monthly runoff (see formula above)
  .[, total_water_without_runoff := 
      (SWE_inst + SoilMoist_RZ_inst + CanopInt_inst)/10] %>%
  .[, total_water := 
      (Qs_month + SWE_inst + SoilMoist_RZ_inst + CanopInt_inst)/10] 
# convert mass to volume: vol = kg / (1000 kg/m3) [Density of water] -> vol = (1/1000) m3
# Therefore: kg/m2 -> (1/1000) [m3/m2] -> (1/10) [cm]


gldas_dates =
  copy(gldas_clsm[,.(time)]) %>% distinct() %>%
  .[,ym := substr(as.Date(time, origin = '2000-01-01'), 1, 7)]  

gldas_clsm = 
  gldas_clsm %>% merge(gldas_dates, by = 'time') %>%
  .[ym %in% grace_dates$ym] 


###########With runoff
gldas_clsm_spread = 
  gldas_clsm[,.(total_water, lon, lat, ym)] %>%
  dcast(lon + lat ~ ym, value.var = "total_water") 

#Estimate the new baseline (2002-2020)
bsl_colnum = 164 #199 = 2020; 165 = 2017; 118 = 2012

#Estimate the new baseline land water content
mean_gldas = apply(gldas_clsm_spread[, 3:(bsl_colnum-1)], 1, mean) 

#Remove the new baseline (currently not using data from GRACE-FO)
gldas_lwc_base = 
  gldas_clsm_spread[, 3:198] - mean_gldas

gldas_clsm_base_rem = 
  cbind(gldas_clsm_spread[,1:2], gldas_lwc_base) %>%
  merge(grace_ref[,.(lon,lat,cell_id)], by = c("lat", "lon"), all.x = F)

#Test the grid by plotting it
test_grid = 
  gldas_clsm_base_rem[,c(2,1,3)] %>% 
  rasterFromXYZ(crs = crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')) 
plot(test_grid)

fwrite(gldas_clsm_base_rem,
       paste0(proj_dir, 'CLSM/GLDAS_CLSM_1degree_BSL2017_SWS_with_runoff_220725.csv'))
