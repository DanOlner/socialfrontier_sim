#Check geographies
#How many OAs in Sheffield?
library(sf)
library(tidyverse)

#These don't have names in.
shefoa <- st_read('C:/Data/MapPolygons/England/2011/England_outputareas_2011_gen/England_oa_2011_gen.shp')

sheflsoa <- st_read('C:/Data/MapPolygons/England/2011/England_lsoa_2011_gen_clipped/england_lsoa_2011_gen_clipped.shp')

sheflsoa$name[grepl(x = sheflsoa$name, pattern = 'Sheffield')]

sheflsoa <- sheflsoa %>% 
  filter(grepl(x = sheflsoa$name, pattern = 'Sheffield'))

plot(sheflsoa)  



rothlsoa <- st_read('C:/Data/MapPolygons/England/2011/England_lsoa_2011_gen_clipped/england_lsoa_2011_gen_clipped.shp')

rothlsoa <- rothlsoa %>% 
  filter(grepl(x = rothlsoa$name, pattern = 'Rotherham'))



