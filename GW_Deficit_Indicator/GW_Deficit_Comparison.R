#Compare GW deficit indicators - downscaled vs non-downscaled

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

############################################################################################
###Load  the paths
proj_dir = "~/Dropbox/WB/GRACE_Ensemble/"
pathOut =  "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"
pathIn = "/Users/tejasvi/Dropbox/WB/GRACE-Deficit/"

############################################################################################

files = list.files(pathIn, '230103', full.names = T)

files.names = 
  list.files(pathIn, '230103', full.names = F) 
files.names = gsub("\\..*","",files.names)

#Load the datasets
ggdi.data = lapply(files, fread)
names(ggdi.data) = files.names

mymerge = function(x,y) {merge.data.table(x,y,all=F,
                                          by.x = c('lat', 'lon'),
                                          by.y = c('Lat', 'Lon'))}

#Merge the datasets
ggdi.merge = 
  Reduce(mymerge, ggdi.data) %>%
  mutate(sumDef = ifelse(Def.bin24_150.x  == 0 & Def.bin24_150.y == 0, 1, 
                         ifelse(Def.bin24_150.x  == 1 & Def.bin24_150.y == 0, 2,
                                ifelse(Def.bin24_150.x  == 0 & Def.bin24_150.y == 1, 3, 4))))


############################################################################################

#Overall
table(ggdi.merge$Def.bin24_150.x,
      ggdi.merge$Def.bin24_150.y,
      dnn = c("Dscl", "non-dscl"), useNA = 'no') %>% prop.table()


plot3 = 
  rasterFromXYZ(ggdi.merge[,.(lon, lat, sumDef)])

plot3 = 
  crop(plot3, extent(-180, 180, -60, 60))  %>%
  rast()

cls <- data.frame(id=1:4, cover=c("No Deficit", "Deficit - Downscaled only",
                                  "Deficit - non-Downscaled only",
                                  "Deficit - Both"))

levels(plot3) <- cls


#color_pal = c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')


#Save the plot for future reference
plot_name =
  paste0("~/Dropbox/WB/GRACE-Deficit/Figures/", "GW_Deficit_24_month_binary_comparison.png")

png(plot_name, width = 1250, height = 500)

plot(plot3, col = c('#ffff30','#beaed4','#fdc086','#7fc97f'),
     plg=list(cex=1.2), mar = c(3.1, 3.1, 2.1, 17.1))

dev.off()

############################################################################################


