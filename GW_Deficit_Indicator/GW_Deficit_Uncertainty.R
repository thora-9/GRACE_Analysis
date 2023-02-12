#Produce groundwater deficit indicator datasets

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

############################################################################################
###Load all the GWS datasets
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"
pathOut =  "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"
pathIn = "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"
pathData = '/Users/tejasvi/Dropbox/gwflagship_GRACEdownscaling/Downscaled TWS_GWS v2/'
pathTyp = '/Users/tejasvi/Dropbox/WB/Typology/'

############################################################################################

#Load the studyregions
med.v = st_read(paste0(pathData, 'studyregions/med.shp'))
sa.v = st_read(paste0(pathData, 'studyregions/sa.shp'))
saf.v = st_read(paste0(pathData, 'studyregions/saf.shp'))

studyRegion = 
  st_union(med.v, sa.v) %>%
  st_union(saf.v)


#######
#Aquifer Typology
## Vector version
aq.v <-
  st_read(paste0(pathTyp, "aqtyp_dissolved.gpkg")) %>%
  # dplyr::select(aqtyp) %>%
  # st_union()
  mutate(aqtyp=factor(aqtyp, levels=c("Major Alluvial","Complex","Karstic","Local/Shallow")))

problem =
  aq.v[84800,]

#######
files = list.files(pathIn, '230128', full.names = T)

files.names = 
  list.files(pathIn, '230128', full.names = F) 
files.names = gsub("\\..*","",files.names)

#Load the datasets
ggdi.data = fread(files)

#Convert the data to sf point
ggdi <- st_as_sf(ggdi.data, coords = c("lon", "lat"), crs = crs(studyRegion))


ds1 = 
  tidync(paste0(pathData, 'GWS_mean_05deg.nc')) %>%
  hyper_tibble() %>% as.data.table() %>%
  dplyr::rename(GWSA_CLSM_En = GWSA_CLSM, GWSA_Noah_En = GWSA_Noah)

ds2 = 
  tidync(paste0(pathData, 'GWS_CSR_05deg.nc')) %>%
  hyper_tibble() %>% as.data.table() %>%
  dplyr::rename(GWSA_CLSM_CSR = GWSA_CLSM, GWSA_Noah_CSR = GWSA_Noah)

ds3 = 
  tidync(paste0(pathData, 'GWS_GSFC_05deg.nc')) %>%
  hyper_tibble() %>% as.data.table() %>%
  dplyr::rename(GWSA_CLSM_GSFC = GWSA_CLSM, GWSA_Noah_GSFC = GWSA_Noah)

 
ds4 = 
  tidync(paste0(pathData, 'GWS_JPL_05deg.nc')) %>%
  hyper_tibble() %>% as.data.table() %>%
  dplyr::rename(GWSA_CLSM_JPL = GWSA_CLSM, GWSA_Noah_JPL = GWSA_Noah)


#Create a date sequence
date.seq = 
  seq(as.Date("2003/2/1"), as.Date("2021/9/1"), "month") %>%
  as.data.table() %>%
  rownames_to_column() %>%
  .[, rowname := as.integer(rowname)]
colnames(date.seq) = c('rowname', 'date')

#Merge the dates using the rowname column
ds1 = 
  ds1 %>%
  merge(date.seq, by.x = 'time', by.y = 'rowname', all.x = T) %>%
  .[, ':='(yearmon = substr(date, 1, 7),
           year = substr(date, 1, 4),
           month = substr(date, 6, 7),
           GWSA_mean_En = (GWSA_CLSM_En+GWSA_Noah_En)/2,
           cell_id = paste0(lat, lon))] %>%
  merge(ds2, by = c('lat', 'lon', 'time'), all.x=T) %>%
  merge(ds3, by = c('lat', 'lon', 'time'), all.x=T) %>%
  merge(ds4, by = c('lat', 'lon', 'time'), all.x=T)
  

gws.unique = 
  ds1 %>% 
  dplyr::select(lat, lon) %>%
  distinct() %>%
  rownames_to_column() 

trends = 
  ds1 %>%
  filter(year<2021) %>%
  merge(gws.unique, by = c('lat', 'lon'), all.x=T) %>%
  group_by(lat, lon, year) %>%
  summarise(GWSA_CLSM_En = mean(GWSA_CLSM_En, na.rm = T),
            GWSA_Noah_En = mean(GWSA_Noah_En, na.rm = T),
            GWSA_mean_En = mean(GWSA_mean_En, na.rm = T),
            GWSA_CLSM_CSR = mean(GWSA_CLSM_CSR, na.rm = T),
            GWSA_Noah_CSR = mean(GWSA_Noah_CSR, na.rm = T),
            GWSA_CLSM_GSFC = mean(GWSA_CLSM_GSFC, na.rm = T),
            GWSA_Noah_GSFC = mean(GWSA_Noah_GSFC, na.rm = T),
            GWSA_CLSM_JPL = mean(GWSA_CLSM_JPL, na.rm = T),
            GWSA_Noah_JPL = mean(GWSA_Noah_JPL, na.rm = T))
  #group_by(lat, lon) %>%
  #group_modify(~as.data.frame(t(zyp.yuepilon(.x$GWS.year)[c(2,6)]))) %>%
  #as.data.table()

df_names = c('GWSA_CLSM_En', 'GWSA_Noah_En', 'GWSA_mean_En', 'GWSA_CLSM_CSR',
             'GWSA_Noah_CSR', 'GWSA_CLSM_GSFC', 'GWSA_Noah_GSFC', 'GWSA_CLSM_JPL',
             'GWSA_Noah_JPL')
i = 1

for (i in 1:length(df_names)){
  print(i)
  cur_var = df_names[i]
  neg_sig = paste0('neg_sig',i)
  tr.name = paste0('trend',i)
  sig.name = paste0('sig',i)
  
  cur_sub = 
    trends[,c('lat','lon', as.character(cur_var))] %>%
    group_by(lat, lon) %>%
    summarise(across(all_of(cur_var), ~as.data.frame(t(zyp.yuepilon(.x)[c(2,6)])))) %>%
    as.data.table() %>%
    rename(!!tr.name := 3, !!sig.name := 4) %>%  #Use !! and := to have strings on the LHS and RHS
    mutate(!!neg_sig := ifelse(.[[3]] < 0 & .[[4]] <= 0.05, 1, 0))
  
  if(i == 1) {
    out.df = cur_sub
  } else {
    out.df = 
      out.df %>%
      merge(cur_sub, by = c('lat', 'lon'), all.x = T)
  }
  
}


out.df2 = 
  out.df %>%
  dplyr::select(lat, lon, contains("neg_")) %>%
  rowwise() %>% 
  dplyr::mutate(total = sum(c_across(starts_with("neg_")), na.rm = T)) %>%
  as.data.table()


cur_rast = 
  out.df2 %>% 
  dplyr::select('lon', 'lat', 'total') %>%
  rasterFromXYZ()

cur_rast[values(cur_rast)==0,] = NA

cur_rast = 
  crop(cur_rast, extent(-180, 180, -60, 60)) %>%
  as.data.frame(xy=T, na.rm = T)


temp.rast = 
  out.df2 %>% 
  dplyr::select('lon', 'lat', 'total') %>%
  rasterFromXYZ()

#Clip the aquifer typlogy raster
aq.v.clipped = 
  aq.v[-84800, ] %>%
  st_make_valid()

aq.v.dscl = 
  st_crop(aq.v.clipped, temp.rast)

############################################################################################

#Merge Trend/Negative 
#Overall
table(ggdi.data$Def.bin24_100,
      ggdi.data$neg_sig,
      dnn = c("Def_10_20", "neg_sig"), useNA = 'no') %>%
  prop.table()

nrow(ggdi.data[neg_sig==1])/nrow(ggdi.data)
nrow(ggdi.data[Def.bin24_100==1])/nrow(ggdi.data)
#

#Aquifer Type
tab.aqtyp = 
  table(ggdi.data$Def.bin24_100,
        ggdi.data$aqtyp_max,
        dnn = c("Def_10_20", "aquifer type"), useNA = 'no') %>% as.data.frame() %>%
  group_by(aquifer.type) %>%
  mutate(total = sum(Freq)) %>%
  mutate(prop = Freq/total) %>%
  filter(Def_10_20 == 1)

tab.aqtyp = 
  table(ggdi.data$neg_sig,
        ggdi.data$aqtyp_max,
        dnn = c("Def_10_20", "aquifer type"), useNA = 'no') %>% as.data.frame() %>%
  group_by(aquifer.type) %>%
  mutate(total = sum(Freq)) %>%
  mutate(prop = Freq/total) %>%
  filter(Def_10_20 == 1)

temp.rast = 
  ggdi.data %>% 
  dplyr::select('lon', 'lat', 'neg_sig') %>%
  rasterFromXYZ()

#Percent negative trends in each TBA
TBA_neg = 
  raster::extract(temp.rast, 
                  trans.all)

len1 = sapply(TBA_neg, function(x){length(x)})
sumNA = sapply(TBA_neg, function(x){sum(is.na(x))})
sumNA.n = sapply(TBA_neg, function(x){sum(!is.na(x))})
sum = sapply(TBA_neg, function(x){sum(x, na.rm = T)})


TBA_neg = 
  cbind(trans.all, len1, sumNA, sumNA.n, sum) %>%
  mutate(out.SR = sumNA/len1) %>%
  filter(out.SR != 1) %>%
  mutate(atleast1 = ifelse(sum>0, 1, 0))

100*sum(TBA_neg$sum>0)/nrow(TBA_neg)

############################################################################################
## SETUP FOR BAR AND BOXPLOTS
scale.aq <- scale_fill_manual(values=c("#44546a", "#70ad47", "#b7ff4b", '#ffc000'), name = "Aquifer type", 
                              guide = guide_legend(order = 1))
theme.blank <-   theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())

plot2 = 
  aq.v.dscl %>%
  drop_na() %>%
  ggplot() +
  geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
  scale.aq +
  new_scale_fill() +
  geom_tile(data = cur_rast, aes(x = x, y = y, fill = total), alpha = 0.8) +
  scale_fill_gradient(low = '#FFD2C1', high = '#A62C2B',
                      breaks = c(1,9),labels = c("Low", "High"), 
                      name = 'Certainty of Depletion') +
  new_scale_fill() +
  geom_sf(data = studyRegion, fill = NA, size = 0.8) +
  new_scale_fill() +
  theme_bw() + ylim(-34.83427, 37.5) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title=element_text(size=16, vjust = 2, face="bold"),
        legend.title.align=0.5,
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathOut, 'Figures/Trend_uncertainty_230208', '.png'), plot=plot2,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')

# ggsave(paste0(pathOut, 'Figures/Trend_uncertainty_230208', '.eps'), plot=plot2,
#        scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm', device=cairo_ps, 
#        fallback_resolution=600)

plot3 = 
  aq.v.dscl %>%
  drop_na() %>%
  ggplot() +
  geom_sf(aes(fill = aqtyp), lwd = 0, alpha = 0.5) + 
  scale.aq +
  new_scale_fill() +
  geom_tile(data = cur_rast, aes(x = x, y = y, fill = total), alpha = 0.8) +
  scale_fill_gradient(low = '#FFD2C1', high = '#A62C2B',
                      breaks = c(1,9),labels = c("Low", "High"), 
                      name = 'Certainty of Depletion') +
  new_scale_fill() +
  geom_sf(data = studyRegion, fill = NA, size = 0.8) +
  new_scale_fill() +
  geom_sf(data = ggdi %>% filter(Def.bin24_100==1), aes(colour = 'Deficit 10/20'), size = 0.2) +
  scale_color_manual(values = 'black', name = NULL, guide = guide_legend(order = 2)) +
  theme_bw() + ylim(-34.83427, 37.5) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.title=element_text(size=16, vjust = 2, face="bold"),
        legend.title.align=0.5,
        #legend.position="right",
        legend.position = c(.15, .26),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathOut, 'Figures/Trend_uncertainty_wdef_230208', '.png'), plot=plot3,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')

# ggsave(paste0(pathOut, 'Figures/Trend_uncertainty_wdef_230208', '.eps'), plot=plot3,
#        scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm', device="eps")

####################################################################
#Test to make sure the Ensemble mean gives the same result
####################################################################

cur_rast = 
  out.df2 %>% 
  dplyr::select('lon', 'lat', 'neg_sig3') %>%
  rasterFromXYZ()

cur_rast[values(cur_rast)==0,] = NA

cur_rast = 
  crop(cur_rast, extent(-180, 180, -60, 60)) %>%
  as.data.frame(xy=T, na.rm = T) %>%
  mutate(deficit=factor(neg_sig3)) %>%
  as.data.table()

levels(cur_rast$deficit) = label


plot2 = 
  ggplot() +
  geom_sf(data = studyRegion, fill = NA, size = 0.8) +
  new_scale_fill() +
  geom_tile(data = cur_rast, aes(x = x, y = y, fill = deficit), alpha = 0.8) +
  scale_fill_manual(values='#d6604d', name = 'Deficit') +
  theme_bw() + ylim(-34.83427, 37.5) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 
