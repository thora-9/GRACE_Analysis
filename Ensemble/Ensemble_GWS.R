#Create final GRACE dataset

###Load all the GWS datasets
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"

gws.path = list.files(paste0(proj_dir, "Outputs/GWS/"), '.csv', full.names = T)

gws_all = lapply(gws.path, fread)

#Get the names for each of them
filePath = fread(paste0(proj_dir, 'FileSummary.csv'))
#Check filepaths in the csv and folder align before adding names
gws.path == filePath[ID %in% c(101:110)]$FilePath

names(gws_all) = filePath[ID %in% c(101:110)]$Name

#Convert each dataset into long format/add new cell_id column

for(i in 1:length(gws_all)){
  cur_df = 
    gws_all[[i]] %>% 
    melt(id.vars = c("lon", "lat"),
         measure.vars = 4:(ncol(.)),
         variable.name = "ym", value.name = names(gws_all[i])) %>%
    mutate(ID = paste0(lon, lat))
  
  if(i == 1) {
    gws_out = cur_df
  } else{
    gws_out = merge(gws_out, 
                    cur_df[,c('ID', 'ym', names(gws_all[i])), with=FALSE], by = c('ID', 'ym') ,all=T)
  }
  
}

################################################################################################
#Use a example gws_wide file to merge spatially with World Bank regions
gws_spatial = 
  gws_all[[1]] %>%
  mutate(ID = paste0(lon, lat)) %>% 
  dplyr::select(lon, lat, ID) %>%
  st_as_sf(coords = c("lon", "lat"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs")

#Add the World Bank regions

####Load World Regions
wb_regions = 
  st_read(paste0(proj_dir,"Spatial Files/WB_Regions/WB_countries_Admin0_10m.shp")) %>%
  dplyr::select(WB_NAME, ISO_A2, ISO_A3, ISO_N3, TYPE) %>%
  filter(TYPE != 'Dependency') %>%
  st_make_valid()

wb_regions_ns = 
  wb_regions %>% as.data.table() %>% dplyr::select(-geometry) %>% distinct()

#Merge country data with GRACE
gws_out_country = 
  gws_spatial %>%
  st_make_valid() %>%
  st_join(wb_regions) %>%
  st_drop_geometry()

#Merge back with the gws_all long data frame
gws.final = 
  gws_out %>%
  merge(gws_out_country, by = 'ID', all = T)

#########################################################################
# fwrite(gws.final,
#        paste0(proj_dir, "Outputs/Ensembles/GRACE_GWS_Ensemble_1degree_220828.csv"))
# 


#############Tests

gws_in = fread(paste0(proj_dir, "Outputs/Ensembles/GRACE_GWS_Ensemble_1degree_220828.csv"))

gws_all[[4]]

View(gws.final[ID == '-178.566.5'])
