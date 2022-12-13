#Downscaled GRACE analysis
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
library(stars)
library(ggnewscale)
library(ggpattern)
library(gridExtra)
library(ggrepel)


proj_dir = "/Users/tejasvi/Dropbox/WB/GRACE_Ensemble/"
pathTyp = '/Users/tejasvi/Dropbox/WB/Typology/'
pathData = '/Users/tejasvi/Dropbox/gwflagship_GRACEdownscaling/Downscaled TWS_GWS v2/'


gws_mean_05 = 
  tidync(paste0(pathData, 'GWS_mean_05deg.nc')) %>%
  hyper_tibble() %>% as.data.table()

#Create a date sequence
date.seq = 
  seq(as.Date("2003/2/1"), as.Date("2021/9/1"), "month") %>%
  as.data.table() %>%
  rownames_to_column() %>%
  .[, rowname := as.integer(rowname)]
colnames(date.seq) = c('rowname', 'date')

#Merge the dates using the rowname column
gws_mean_05 = 
  gws_mean_05 %>%
  merge(date.seq, by.x = 'time', by.y = 'rowname', all.x = T) %>%
  .[, ':='(yearmon = substr(date, 1, 7),
           GWSA_mean = (GWSA_CLSM+GWSA_Noah)/2)] 

#Convert to wide
gws_05_wide = 
  gws_mean_05 %>% 
  dplyr::select(GWSA_mean, lat, lon, date) %>%
  tidyr::spread(key = date, value = GWSA_mean) %>%
  mutate(cell_id = paste0(lon, lat))

#Convert to spatial object
gws_05.v = 
  st_as_sf(gws_05_wide, coords = c("lon", "lat"), 
           crs = "+proj=longlat +datum=WGS84 +no_defs")

#Load the studyregions
med.v = st_read(paste0(pathData, 'studyregions/med.shp'))

sa.v = st_read(paste0(pathData, 'studyregions/sa.shp'))

saf.v = st_read(paste0(pathData, 'studyregions/saf.shp'))

regions_all = 
  bind_rows(med.v, sa.v, saf.v)

continent = 
  st_read(paste0(pathTyp, 'IGRAC_Transboundary/World_Continents.shp')) %>%
  st_transform(crs(regions_all))


## transboundary aquifers
trans.sub = c('AF058','AF064','AF056', 'AF001','AF054','AF063',
              'AS126','AS150','AS079','AS080','AS089','AS090','AS091',
              'AS150','S021','S015','C007')

trans <-
  st_read(paste0(pathTyp, 'IGRAC_Transboundary/2021_TBA_GGIS_utf8_VALID.shp')) %>%
  filter(CODE_2021 %in% trans.sub) %>%
  st_transform(crs(regions_all)) %>%
  st_make_valid() %>%
  distinct()

#-------------------------------------------------------------------------------
# GRACE downscaled map
#-------------------------------------------------------------------------------


#Global Map

plot2 = 
  continent %>%
  drop_na() %>%
  ggplot() +
  geom_sf(aes(geometry = geometry), fill=NA) + 
  geom_sf(regions_all,
          mapping = aes(fill= 'Downscaled GRACE\nData Available')) +
  scale_fill_manual(name = NA, values = 'darkolivegreen2') +
  new_scale_fill() +
  geom_sf_pattern(data = trans, 
                  aes(pattern = 'Transboundary\nAquifer'), colour = 'firebrick4',
                  fill = NA, pattern_colour  = 'black', pattern_alpha = 0.7, 
                  pattern_density = 0.6, pattern_fill = NA,
                  pattern_spacing = 0.005) +
  geom_label_repel(data = trans,
                   aes(label = CODE_2021, geometry = geometry),
                   stat = "sf_coordinates",
                   min.segment.length = 0) +
  scale_pattern_manual(name = NA, values = 'circle') +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        #legend.position="right",
        legend.position = c(.15, .3),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.text=element_text(size=16), legend.title=element_blank(), legend.key.size = unit(0.8, "cm"),
        panel.background = element_blank(), axis.line = element_blank()) 

ggsave(paste0(pathTyp, 'plots/trans_GRACE_Global', '.png'), plot=plot2,
       scale=1.5, dpi=300,width =34.85,height = 18, units = 'cm')


#-------------------------------------------------------------------------------
# GRACE Trends
#-------------------------------------------------------------------------------


#First, convert values from monthly to year
gws_annual = 
  gws_05_wide %>% 
  data.table::melt(id.vars = c("lon", "lat", "cell_id"),
                   measure.vars = 3:214,
                   variable.name = "ymd", value.name = "GWSA_mean") %>% 
  .[, year := substr(ymd, 1,4)] %>%
  .[, an_mean := mean(GWSA_mean), .(lat, lon, year)] %>% 
  .[year<2021] %>%
  .[,.(lat, lon, cell_id, year, an_mean)] %>% distinct() %>%
  dcast(lon + lat + cell_id ~ year, value.var = "an_mean") %>% .[order(cell_id)]

gws_trends = 
  apply(gws_annual[, 4:ncol(gws_annual)], 1, zyp.yuepilon)

gws_trends_estimates = 
  gws_annual[,1:3] %>%
  mutate(trends = gws_trends[2,], sig = gws_trends[6,]) %>%
  mutate(lon = ifelse(lon>180, lon-360, lon))

gws_trends_estimates = 
  gws_trends_estimates %>%
  mutate(trends = ifelse(trends < -2, -2, trends)) %>%
  mutate(trends = ifelse(trends > 2, 2, trends)) %>%
  mutate(cuts = cut(trends,
                    breaks = ()))
  

gws.raster = rasterFromXYZ(gws_trends_estimates[,.(lon, lat, trends)])

color_pal = c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')

ggplot() +
  geom_raster(gws_trends_estimates,
              mapping = aes(x = lon, y = lat, fill = trends)) + 
  scale_fill_gradientn(colours = color_pal) + 
  coord_quickmap()


plot(gws.raster,
     breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2),
     col = color_pal, 
     legend.args = list(text = 'Trend (cm/yr)', side = 4, 
                        font = 2, line = 2.5, cex = 0.8))

