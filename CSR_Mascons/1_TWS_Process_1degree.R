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

plot_trends <- function(R, horiz){
  
  #R = wb_region_list[2]
  
  region = 
    wb_regions %>% 
    filter(WB_REGION == R | FORMAL_EN == R) %>% 
    st_transform(x, crs = st_crs(3857))
  
  if(nrow(region) > 0){
    
    bounds = st_bbox(region) 
    # bounds[bounds > 180] = 179
    # bounds[bounds < -179] = 0
    
    region_trends_temp = 
      gws_trends_spatial %>% st_transform(x, crs = st_crs(3857)) %>%
      st_within(st_make_valid(region)) 
    
    sel_logical = lengths(region_trends_temp) > 0
    
    region_trends = 
      gws_trends_spatial[sel_logical, ] %>% st_transform(x, crs = st_crs(3857)) %>%
      st_transform(x, crs = st_crs('+proj=longlat +datum=WGS84 +no_defs')) %>%
      st_cast("MULTIPOINT") %>%  
      filter(!is.na(trends))# %>% st_shift_longitude()
    
    #Create a grid for interpolation
    region_out = region %>% st_transform(x, crs = st_crs(4326))
    
    bbox <- st_bbox(region_out) #%>% st_transform(x, crs = st_crs(4326)))
    
    grd_template <- expand.grid(
      X = seq(from = bbox["xmin"], to = bbox["xmax"], by = 0.2),
      Y = seq(from = bbox["ymin"], to = bbox["ymax"], by = 0.2) # 20 m resolution
    )
    
    #Projection for Raster
    crs_raster_format <- "+proj=longlat +datum=WGS84 +no_defs"
    
    #Create raster template
    grd_template_raster <- grd_template %>% 
      dplyr::mutate(Z = 0) %>% 
      raster::rasterFromXYZ( 
        crs = crs_raster_format)
    
    # Inverse Distance Weighting
    fit_IDW <- gstat::gstat( # The setup here is quite similar to NN
      formula = trends ~ 1,
      data = as(region_trends, "Spatial"),
      nmax = 10, nmin = 3,
      set = list(idp = 0.5) # inverse distance power
    )
    
    interp_IDW <- interpolate(grd_template_raster, fit_IDW)
    
    ## crop and mask
    r2 <- crop(interp_IDW, extent(region_out))
    r3 <- terra::mask(r2, region %>% st_transform(x, crs = st_crs(4326)))
    
    r3_df = 
      as.data.frame(r3, xy = T) %>% rename(values = var1.pred) %>% filter(!is.na(values)) %>%
      mutate(trends_cut = cut(values, 
                              breaks = c(-Inf, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, Inf))) %>%
      dplyr::select(-trends_cut) %>%
      mutate(values = ifelse(values < -2, -2, values)) %>%
      mutate(values = ifelse(values > 2, 2, values)) %>%
      rasterFromXYZ()
    
    
    color_pal = c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')
    
    out_file = paste0("Plots/GWS_",R,".png")
    out_title = paste('Groundwater Storage Trends in', R,"(2002-2017)")
    
    # test = 
    #   ggplot() + 
    #   geom_sf(data = region_out, col = 'black') +
    #   geom_raster(r3_df, mapping = aes_string(x = 'x', y = 'y', fill = 'values')) +
    #   theme_bw() +
    #   coord_sf(datum = st_crs(region_out)) +
    #   
    
    
    png(file= out_file,
        width=1000, height=800)
    
    #pe <- as.polygons(ext(r3_df), crs = 'WGS84')
    
    plot(r3_df, legend=FALSE,
         breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
         col = color_pal) 
    
    title(main = out_title, line = 1, adj = 0.45, cex.main=2)
    
    
    if(horiz == T) {
      loc.1 = c(0.07, 0.34, 0.2, 0.24) #Change these values to adjust the legend position
    } else if (horiz == F) {
      loc.1 = c(0.2, 0.24, 0.18, 0.45)
    }
    
    if(R %like% 'Yemen'){
      loc.1 = c(0.27, 0.54, 0.2, 0.24) #Change these values to adjust the legend position
    }
    
    plot(r3_df, legend.only = T, smallplot = loc.1, horizontal = horiz,
         breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), lab.breaks = c('< -2', -1.5, -1, -0.5, 0, 0.5, 1, 1.5, ">2"),
         col = color_pal, legend.args = list(text = 'Groundwater Storage\n Trends (cm/yr)', side = 3, 
                                             font = 2, line = 1, cex = 1.5)) 
    
    region = st_transform(region, crs = st_crs(4326))
    plot(region$geometry, add = T, lwd = 1)
    
    dev.off()
  }
}

#Load the datasets

#The GRACE Mascon dataset
tws_mascon =
  tidync(paste0(proj_dir,'GRACE_data/CSR_Mascons/CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02_1degree.nc')) %>%
  activate("D2,D3,D0") %>%
  hyper_tibble() %>% as.data.table() %>%
  .[order(lat, lon)] %>% #ensures the order is the same
  group_by(lat,lon) %>%
  mutate(cell_id = cur_group_id())  %>% as.data.table()


dates =
  tidync(paste0(proj_dir,'GRACE_data/CSR_Mascons/CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc')) %>%
  activate("D1") %>%
  hyper_tibble() %>% as.data.table() 

#2 of the dates don't align with the JPL solutions by 1 day $ 4 days
#Manually add those days to make them align
dates[111,]$time = dates[111,]$time + 1
dates[112,]$time = dates[112,]$time + 4

dates = 
  dates %>%
  mutate(dates = as.Date(time, origin = '2002-01-01')) %>%
  mutate(ym = substr(ymd(dates), 1, 7)) 

dates1 = dates$dates 

#This dates merge works because the 'distinct' command has not been run on
#the dates

########################
######### Analysis
########################
#Add dates to the TWS df; also scale the values 
tws_mascon = 
  tws_mascon %>% group_by(cell_id) %>% 
  mutate(ym = dates1) %>% as.data.table()

#Convert from long to wide format
#Two of the months have 2 values that were taken (2011-10 and 2015-04) 
tws_lwe = 
  tws_mascon[,.(lwe_thickness, lon, lat, cell_id, ym)] %>%
  dcast(lon + lat + cell_id ~ ym, value.var = "lwe_thickness") 

#Drop the two ym months
tws_lwe = 
  tws_lwe %>%
  dplyr::select(-148, -116)

colnames(tws_lwe)[4:210] = substr(colnames(tws_lwe)[4:210], 1, 7)


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
       paste0(proj_dir, "GRACE_Data/CSR_Mascons/GRACE_CSR_TWS_1degree_220810.csv"))

