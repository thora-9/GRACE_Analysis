#Transform datasets to 0.5 degree

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


###Load all the GWS datasets
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"

#Load the non-0.5 degree datasets
ensemble = 
  fread(paste0(proj_dir, 
               'Outputs/GWS/9003_GRACE_GWS_TWS_3Mascons_SWS_Ensemble_2002_2020_BSL2017.csv')) %>%
  #mutate(lon = lon + 0.001, lat = lat + 0.001) %>% #Jitter so points fall within a single polygon
  st_as_sf(coords = c("lon", "lat"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_buffer(dist = 0.1)

ensemble.r = 
  fread(paste0(proj_dir, 
               'Outputs/GWS/9003_GRACE_GWS_TWS_3Mascons_SWS_Ensemble_2002_2020_BSL2017.csv')) 
ensemble.r =
  rasterFromXYZ(ensemble.r[,.(lon, lat, `2021-01`)])

#Load the fishnet
fishnet = 
  st_read('/Users/tejasvi/Dropbox/WB/Fishnet_halfdegree/global_fishnet.shp')

#GRACE-ensemble to fishnet
fishnet.GRACE = 
  fishnet %>%
  filter(!row_number() %in% c(18901, 23763, 26074, 42334, 31737)) %>% #Invalid polygon
  st_make_valid() %>%
  st_join(ensemble)

test =
  rasterFromXYZ(fishnet.GRACE[,c('Lon', 'Lat', "2019-12")] %>% st_drop_geometry())

#Get point form data
fishnet.centroid = 
  fishnet.GRACE %>%
  st_centroid()


####Load World Regions
wb_regions = 
  st_read(paste0(proj_dir,"Spatial Files/WB_Regions/WB_countries_Admin0_10m.shp")) %>%
  dplyr::select(WB_NAME, ISO_A2, ISO_A3, ISO_N3, TYPE, WB_REGION) %>%
  filter(TYPE != 'Dependency') %>%
  st_make_valid()

#Merge country data with GRACE
fishnet.centroid = 
  fishnet.centroid %>%
  st_make_valid() %>%
  st_join(wb_regions) 

fishnet.df = 
  fishnet.centroid %>% 
  drop_na() %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  as.data.table()

fishnet.df.long = 
  fishnet.df %>%
  pivot_longer(
       cols = `2002-04`:`2020-12`,
       names_to = "ym", values_to = 'GWS_ensemble') %>%
  dplyr::select(-`2021-01`, -OBJECTID) %>%
  as.data.table()

test =
  rasterFromXYZ(fishnet.df[,.(Lon, Lat, `2021-01`)])

#Write output (spatial)
st_write(fishnet.centroid, paste0(proj_dir, 
             'Outputs/GWS/half-degree/GRACE_GWS_3Mascons_2002_2020_BSL2017_05degree.shp')) 

#Write output (non-spatial)

#Global
fwrite(fishnet.df.long, paste0(proj_dir, 
                          'Outputs/GWS/half-degree/GRACE_GWS_3Mascons_2002_2020_BSL2017_05degree.csv'))

#South-Asia
fwrite(fishnet.df.long[WB_REGION == 'SOA'], paste0(proj_dir, 
                          'Outputs/GWS/half-degree/GRACE_GWS_4Ensemble_2002_2020_BSL2017_05degree_SA.csv'))


#A visual comparison between fishnet, my GRACE grid and Di Long's grid suggests that 
#some of the coastal pixels are lost. The loss seems to inconsistent between the two GRACE
#products, and overall, the loss is greater in the downscaled product.


test = 
  fread('/Users/tejasvi/Dropbox/WB/GRACE_Analysis/GLDAS_2002_2017_SWS_with_runoff.csv')

test2 = 
  fread('/Users/tejasvi/Dropbox/WB/GRACE_Analysis/GWS_2002_2017_with_runoff.csv')

sum(is.na(test$`2003-02`))
sum(is.na(test2$`2003-02`))

sum(is.na(fishnet.centroid$`2002-04`))

sum(test$lat %in% fishnet.centroid$Lat)
sum(test$lon %in% fishnet.centroid$Lon)

fwrite(test, paste0(proj_dir, 'test.csv'))
