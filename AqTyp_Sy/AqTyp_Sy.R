#Link aquifer types in India with specific yield from CGWB

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
library(stars)
library(fasterize)
library(foreign)
library(RcppRoll)
library(tmap)
library(ggnewscale)
library(readr)
library(strex)

############################################################################################
###Load all the GWS datasets
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"
pathOut =  "/Users/tejasvi/Dropbox/WB/India_Aquifer_Sy/"
pathIn = "/Users/tejasvi/Dropbox/WB/India_Aquifer_Sy/"
############################################################################################

aq.ty =
  st_read(paste0(pathIn, 'aquifers_all.shp')) 

aq.df = 
  aq.ty %>% st_drop_geometry %>%
  as.data.table() 

sy.extr = stringr::str_extract_all(aq.df$YEILD__, '\\d+([.,]\\d+)?') 

#Get the max values stored in the Yield field
mx <- max(lengths(sy.extr))
out <- do.call(rbind.data.frame, lapply(sy.extr, `length<-`, mx))
names(out) = c('n1', 'n2')

aq.df = 
  aq.df %>%
  mutate(sy1 = out$n1,
         sy2 = out$n2) %>%
  mutate(sy1 = as.numeric(sy1), sy2 = as.numeric(sy2)) %>%
  rowwise() %>%
  mutate(sy_final = mean(c(sy1, sy2), na.rm = T)) %>%
  group_by(NEWCODE43) %>% #Group by aquifer group code
  mutate(sy_group = round(mean(sy_final, na.rm = T))) %>% #Fill the aquifers missing values with group value
  mutate(sy_final = ifelse(sy_final %in% c(NaN, NA, ""), sy_group, sy_final)) %>%
  group_by(AQUIFER) %>%
  mutate(sy_aq = round(mean(sy_final, na.rm = T))) %>% #Fill the aquifers missing values with group value
  mutate(sy_final = ifelse(sy_final %in% c(NaN, NA, ""), sy_aq, sy_final))
  
#Method - 
#a) Use the max value stated when sy is stored as 'upto x%'
#b) When a numerical range is provided, take the mean
#c) Fill with group mean, when a particular aq_type is missing a value
#d) Fill with aquifer mean, when both group and individual values are missing

aq.df.cat = 
  aq.df %>%
  dplyr::select(OBJECTID, aq_code = NEWCODE43, aq = AQUIFER, state = STATE, sy_OG = YEILD__, sy1, sy2, sy_final)

#Merge back to the shape file
aq.sy.final = 
  aq.ty %>% dplyr::select(OBJECTID) %>%
  merge(aq.df.cat, by = 'OBJECTID', all.x = T) %>%
  dplyr::rename(id = OBJECTID) %>%
  st_transform(4326)


#Seems like the original aquifer data is missing some Sandstone entries 
#Work around create a raster and fill the missing polygons

#Create a rasterized version of the layer
base.r = raster() %>% crop(aq.sy.final)
res(base.r) = 0.01 

#Fasterize the aq.sy
aq.sy.r = fasterize(aq.sy.final, base.r, 'sy_final')


############################################################################################
#Write Output
#st_write(aq.sy.final, paste0(pathOut, 'IN_aq_sy.shp'), append = F)
#fwrite(aq.df.cat, paste0(pathOut, 'aq_df_cat.csv'))  
#writeRaster(aq.sy.r, paste0(pathOut, 'IN_aq_sy.tif'), overwrite=TRUE)
############################################################################################

scale.aq <- scale_fill_manual(values=c("#44546a", "#70ad47", "#b7ff4b", '#ffc000'), name = "Aquifer type")
theme.blank <-   theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())


plot1 = 
  aq.sy.final %>%
  ggplot() +
  geom_sf(aes(fill = sy_final), lwd = 0, alpha = 1) + 
  scale_fill_viridis_b(name = "Specific Yield (%)\n", option = "D",
                       breaks = c(0,2,4,7,10,13,15)) +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.75, .2),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title = element_text(size=18),
        legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathOut, 'IN_aq_sy', '.png'), plot=plot1,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')



