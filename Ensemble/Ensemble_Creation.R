#Ensemble Creation
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

###Create Groundwater Storage Anomaly
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"

filePath = fread(paste0(proj_dir, 'FileSummary.csv'))

#Select the GLDAS Solution (Serves to subset the GRACE datasets)
sws_id = 22
validGridPoints = 
  fread(filePath[Type=='GLDAS' & ID == sws_id]$FilePath) %>%
  .[order(lon,lat)] %>%
  mutate(ID2 = paste0(lon, lat))

#Select the GRACE Solutions - JPL
grace_id = 11
tws_anomaly1 = 
  fread(filePath[Type=='GRACE' & ID == grace_id]$FilePath) %>%
  mutate(ID2 = paste0(lon, lat)) %>%
  dplyr::filter(ID2 %in% validGridPoints$ID2) %>% 
  .[order(lon,lat)] %>%
  melt(id.vars = c("lon", "lat", "ID2"),
       measure.vars = 4:(ncol(.)-1),
       variable.name = "ym", value.name = "tws") 


grace_id = 12
tws_anomaly2 = 
  fread(filePath[Type=='GRACE' & ID == grace_id]$FilePath) %>%
  mutate(ID2 = paste0(lon, lat)) %>%
  dplyr::filter(ID2 %in% validGridPoints$ID2) %>% 
  .[order(lon,lat)] %>%
  melt(id.vars = c("lon", "lat", "ID2"),
       measure.vars = 4:(ncol(.)-1),
       variable.name = "ym", value.name = "tws") 

grace_id = 13
tws_anomaly3 = 
  fread(filePath[Type=='GRACE' & ID == grace_id]$FilePath) %>%
  mutate(ID2 = paste0(lon, lat)) %>%
  dplyr::filter(ID2 %in% validGridPoints$ID2) %>% 
  .[order(lon,lat)] %>%
  melt(id.vars = c("lon", "lat", "ID2"),
       measure.vars = 4:(ncol(.)-1),
       variable.name = "ym", value.name = "tws") 
  

tws_anomaly =
  merge(tws_anomaly1, tws_anomaly2[,.(ID2, ym, tws2 = tws)], all=TRUE, by = c('ID2', 'ym')) %>%
  merge(tws_anomaly3[,.(ID2, ym, tws3 = tws)], all=TRUE, by = c('ID2', 'ym')) %>%
  dplyr::mutate(tws_ensemble = (tws+tws2+tws3)/3)

tws_anomaly_wide = 
  tws_anomaly[,.(lon, lat, cell_id = ID2, ym, tws_ensemble)] %>%
  spread(ym, tws_ensemble)

fwrite(tws_anomaly_wide,
       paste0(proj_dir, "Outputs/Ensembles/GRACE_TWS_Ensemble_1degree_220828.csv"))


#Load the sws anomalies - 1
sws_id = 21
sws_anomaly1 = 
  fread(filePath[Type=='GLDAS' & ID == sws_id]$FilePath) %>%
  .[order(lon,lat)] %>%
  mutate(ID2 = paste0(lon, lat)) %>%
  melt(id.vars = c("lon", "lat", "ID2"),
       measure.vars = 3:(ncol(.)-2),
       variable.name = "ym", value.name = "sws") 


#Load the sws anomalies - 2
sws_id = 22
sws_anomaly2 = 
  fread(filePath[Type=='GLDAS' & ID == sws_id]$FilePath) %>%
  .[order(lon,lat)] %>%
  mutate(ID2 = paste0(lon, lat)) %>%
  melt(id.vars = c("lon", "lat", "ID2"),
       measure.vars = 3:(ncol(.)-2),
       variable.name = "ym", value.name = "sws") 

#Load the sws anomalies - 3
sws_id = 23
sws_anomaly3 = 
  fread(filePath[Type=='GLDAS' & ID == sws_id]$FilePath) %>%
  .[order(lon,lat)] %>%
  mutate(ID2 = paste0(lon, lat)) %>%
  melt(id.vars = c("lon", "lat", "ID2"),
       measure.vars = 3:(ncol(.)-2),
       variable.name = "ym", value.name = "sws") 

sws_anomaly =
  merge(sws_anomaly1, sws_anomaly2[,.(ID2, ym, sws2 = sws)], all=TRUE, by = c('ID2', 'ym')) %>%
  merge(sws_anomaly3[,.(ID2, ym, sws3 = sws)], all=TRUE, by = c('ID2', 'ym')) %>%
  dplyr::mutate(sws_ensemble = (sws+sws2+sws3)/3)


sws_anomaly_wide = 
  sws_anomaly[,.(lon, lat, ym, sws_ensemble, cell_id = ID2)] %>%
  spread(ym, sws_ensemble) %>%
  dplyr::select(lon, lat, 4:199, cell_id) #To ensure same format as other GLDAS tables

fwrite(sws_anomaly_wide,
       paste0(proj_dir, "Outputs/Ensembles/GRACE_SWS_Ensemble_1degree_220828.csv"))

##############################################
#Create the gws anomaly here
gws_anomaly = 
  tws_anomaly[,.(ID2,lon,lat,ym,tws_ensemble)] %>%
  merge(sws_anomaly[,.(ID2,ym,sws_ensemble)], by = c('ID2', 'ym')) %>%
  mutate(gws_ensemble = tws_ensemble - sws_ensemble)

gws_anomaly_wide = 
  gws_anomaly[,.(lon,lat, cell_id = ID2, ym,gws_ensemble)] %>%
  spread(ym, gws_ensemble)


