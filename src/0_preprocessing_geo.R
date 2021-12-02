#
# Compute pairwise distances between stations
# identified by lon-lat coordinates
#
# Distances are computed as minimal distances
# between two locations without passing by lands
#
# Jean Olivier Irisson, modified by Ophélie Da Silva
# December, 2021
#-------------------------------------------------- 

# load libraries
library("tidyverse")
library("gdistance")
library("ggrepel")
library("sp")


# clean environment
rm(list=ls())
graphics.off()

source("src/Rsources.R")

# read data
st <- read.table("data/raw/coord.txt", 
                 header = T, stringsAsFactors = F) %>%
  rename(id=Station, lon=Longitude, lat=Latitude) %>% 
  mutate(id = gsub("^TARA_([0-9]+)$", "\\1", id) %>% as.numeric) %>% 
  dplyr::select(id, lon, lat)

coast <- read_csv("data/raw/coastline_medit.csv", col_types=cols())

#-----------------------------------------------------
# Main
#-----------------------------------------------------

# convert coastline to SpatialPolygons
cl <- split(coast, cumsum(is.na(coast$lon)))
cl <- lapply(cl, na.omit)
cl <- cl[-length(cl)]
for (i in 1:length(cl)) {
  cl[[i]] <- Polygons(list(Polygon(cl[[i]])), ID=i)
}

crs <- CRS("+proj=longlat +datum=WGS84")
cp <- SpatialPolygons(cl, proj4string=crs)

# rasterise it with reasonnable precision
res <- 8
cr <- raster(extent(cp), nrow=70*res, ncol=230*res, crs=crs)
cr <- rasterize(cp, cr, background=0)
cr[cr>=1] <- 1
cr <- (1 - cr) # invert it

# compute conductance matrix------------------------------------------------------------
# the conductance function f(x) defines the conductance from x[1] to x[2]
# here = 1 when arriving in sea, 0 when arriving on land = the value of x[2]
landing <- function(x) { x[2] }
cond <- transition(cr*100, transitionFunction = landing, 8)

# and divide conductance by actual geographic distance between cell centers, to get minimum spanning distances
cond <- geoCorrection(cond, type="c", multpl=FALSE, scl=FALSE)
#---------------------------------------------------------------------------------------

## Compute pairwise distances among stations
# get coordinates of populations
coords <- dplyr::select(st, -id) %>% as.matrix()

# compute distances
geo_dist <- costDistance(cond, coords) / 10 # /1000 -> to convert in km

# convert them to long dataframe
as.dist.df <- function(df, labs = st$id) {
  df %>% 
    as.matrix() %>% 
    `colnames<-`(labs) %>% 
    `row.names<-`(labs) %>% 
    as.data.frame() %>% 
    mutate(from = row.names(.)) %>% 
    gather(key = to, value = dist, -from) %>% 
    return()
}

geo_df <- as.dist.df(geo_dist) %>% rename(geo = dist)

# Formating station names in environmental table
geo_df %>% 
  mutate(from = formating_name(from),
         to = formating_name(to)) -> geo_df

# write data to disk
write.table(geo_df, file="data/out/geo-distances.txt", 
            sep="\t", row.names=F, quote=F)