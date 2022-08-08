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

#Select the GLDAS Solution
sws_id = 21
sws_anomaly = 
  fread(filePath[Type=='GLDAS' & ID == sws_id]$FilePath) %>%
  .[order(lon,lat)]

#Select the GRACE Solution
grace_id = 11
tws_anomaly = 
  fread(filePath[Type=='GRACE' & ID == grace_id]$FilePath) %>%
  dplyr::filter(cell_id %in% sws_anomaly$cell_id) %>% 
  .[order(lon,lat)]

combo_name = paste(filePath[ID==grace_id]$Name, filePath[ID==sws_id]$Name, sep='_')

#Checks to ensure the datasets align 
identical(sws_anomaly$cell_id, tws_anomaly$cell_id)
identical(colnames(sws_anomaly)[3:193], colnames(tws_anomaly)[4:194]) #Subset till 2020-Dec
identical(dim(tws_anomaly[ ,4:194]), dim(sws_anomaly[ ,3:193]))

gws_anomaly = 
  tws_anomaly[ ,4:194] - sws_anomaly[ ,3:193]

gws_anomaly = cbind(tws_anomaly[,1:3], gws_anomaly)

#WRITE OUTPUT CSV
fwrite(gws_anomaly, 
       paste0(proj_dir, "Outputs/GWS/GRACE_GWS_", combo_name,
              "_2002_2020_BSL2017.csv"))


####Trend Analysis 
#### (also use this as a visual test to see if expected hotspots pop-out)

#Convert monthly to annual
gws_annual = 
  gws_anomaly %>% melt(id.vars = c("lon", "lat", "cell_id"),
                       measure.vars = 4:ncol(gws_anomaly),
                       variable.name = "ym", value.name = "gws") %>% 
  .[, year := substr(ym, 1,4)] %>%
  .[, an_mean := mean(gws), .(lat, lon, year)] %>% 
  .[,.(lat, lon, cell_id, year, an_mean)] %>% distinct() %>%
  dcast(lon + lat + cell_id ~ year, value.var = "an_mean") %>% .[order(cell_id)]

gws_trends = 
  apply(gws_annual[, 4:ncol(gws_annual)], 1, zyp.yuepilon)

gws_trends_estimates = 
  gws_annual[,1:3] %>%
  mutate(trends = gws_trends[2,], sig = gws_trends[6,]) %>%
  mutate(lon = ifelse(lon>180, lon-360, lon))

gws.raster = rasterFromXYZ(gws_trends_estimates[,.(lon, lat, trends)])

#Save the plot for future reference
plot_name =
  paste0(proj_dir, "Outputs/GWS/GRACE_GWS_", combo_name,
         "_2002_2020_BSL2017.png")

png(plot_name, width = 1250, height = 650)

plot(gws.raster, zlim=c(-10,10))

dev.off()

#fwrite(gws_trends_estimates, "Output/GRACE_GWS_annual_trends_220203.csv")
