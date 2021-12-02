#
#
#
# Extracting environmental data of from Medar-Medatlas
#
#
#
# Ophélie DA SILVA
# December, 2021
#----------------------------------------------------


# General setup
options(stringsAsFactors = F)

# Libraries
library("tidyverse")
library("RNetCDF")

# Cleaning
rm(list=ls())
graphics.off()

# Functions
source("src/Rsources.R") # list2df

# Initial variables
env.f = "data/raw/coord.txt"
input = "data/raw/MEDAR-MEDATLAS"
output = "data/out"
var_info.f = paste(output,"medar_var_info.txt",sep="/")

# MAIN -----------------------------------
# Coordinates of the Tara Oceans stations
tara_stations <- read.table(env.f, header = T)
rm(env.f)

# Rounding depth on same intervals than medar/medatlas
tara_stations %>% 
  rename(latitude = Latitude,
         longitude = Longitude,
         depth = Depth.nominal) %>% 
  mutate(depth = round(depth/5)*5) -> tara_stations
depth <- unique(tara_stations$depth)

# Medar/Medatlas data
medar.f = list.files(input) %>% 
  .[grep("clim", .)] %>% 
  paste(input, ., sep="/"); rm(input)

# Coordonnées géographiques les plus proches de la stations TARA
dat.l = list()
sink(var_info.f)

for (f in medar.f) {
  
  ## Reading NetCDF
  system(paste("gunzip -k", f))
  tmp = gsub(".gz$","",f)
  
  nc = open.nc(tmp)
  print.nc(nc)
  dat = read.nc(nc)
  system(paste("rm",tmp))
  close.nc(nc)
  
  rm(nc)
  
  ## Formating data
  dat.l[[f]] = list(depth = dat$depth,
                    lat = dat$lat,
                    lon = dat$lon,
                    var = dat[[4]]) %>%
    list2df %>%
    `names<-`(c(names(.)[1:3], gsub(".*\\.(.+)\\.nc", "\\1", tmp))) %>%
    .[.[,"depth"] %in% depth ,]
  
  rm(dat, tmp)
}
sink()

Reduce(function(x, y) 
  merge(x, y, all=TRUE), dat.l) -> dat.df

rm(f, medar.f, depth, dat.l)

# Extracting value from the nearest grid cell
full_join(dat.df, tara_stations) %>%
  mutate(dist = 
           6371 * acos(sin(lat*pi/180) * sin(latitude*pi/180) + 
                         cos(lat*pi/180)*cos(latitude*pi/180)*cos(lon*pi/180-longitude*pi/180))) -> tmp

aggregate(tmp$dist, by = list(tmp$latitude, tmp$longitude, tmp$depth), min) %>%
  `colnames<-`(c("latitude", "longitude", "depth", "dist")) %>%
  merge(., tmp, all.x = T) %>%
  merge(tara_stations, ., all.x = T) %>%
  select(-lat,-lon, -dist, -depth) %>% 
  `row.names<-`(gsub("^TARA_0{1,2}([0-9]+)$","\\1", .$Station)) %>% 
  rename(station = Station) -> env

rm(tmp, dat.df, tara_stations, list2df)

# Write data
write.table(env, paste(output,"medar_env_tara_med.txt",sep="/"), row.names = F, quote = F)
