
#---------------------------
# Graphics parameters
#---------------------------


# for environmental clusters
color_scale <- c("#4674a0","#d1495b","#edae49","#90bba4","#a3a3a3")

# for geographic bassins and seas
geo_shape <- c(3, 15, 8, 17, 18, 19)


#---------------------------
# Map background
#---------------------------
maptheme <- theme(panel.grid = element_blank()) +
  theme(axis.text = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = "bottom") +
  theme(panel.background = element_rect(fill = "#596673")) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0), 'cm'))
country_shapes = geom_polygon(aes(x = long, y = lat, group = group),
                              data = map_data('world'),
                              fill = "#969799", color = "#515151",
                              size = 0.15)
mapcoords <- coord_fixed(xlim = c(-5, 35), ylim = c(30.5, 44.5))

ggplot() +
  maptheme +
  country_shapes +
  mapcoords -> map_bck

#------------------------------------------------------
# Stations and basin names correspondances
#------------------------------------------------------
bassin <- data.frame(station = c(5,6,7,9,16,18,20,22,23,24,25,26,30) %>% 
                       as.character(),
                       geo_clust = c(rep("Alboran Sea",2),
                                     rep("Western Mediterranean",2),
                                     rep("Ionian Sea",2), 
                                     rep("Tunisian Plateau/Gulf of Sicila",1), 
                                     rep("Ionian Sea",1), 
                                     rep("Adriatic Sea",2), 
                                     rep("Ionian Sea",2), 
                                     "Levantine Sea"))

#------------------------------------------------------------
# For figure 4
#------------------------------------------------------------

df_to_matrix <- function(x) {
  # create distance matrix
  x %>% 
    select(-Organism) %>%
    spread(to,Fst) %>% 
    as.data.frame() %>% 
    `row.names<-`(.$from) %>% 
    select(-from) -> mat
  
  # Add missing row(s)/column(s) to obtain a squared matrix
  mat[,setdiff(rownames(mat), colnames(mat))]=NA
  mat[setdiff(colnames(mat), rownames(mat)),]=NA
  
  # Order rows  
  row.names(mat) %>%
    as.numeric() %>% 
    order()  %>% 
    mat[.,] -> mat
  
  # Order columns
  colnames(mat) %>%
    as.numeric() %>% 
    order()  %>% 
    mat[,.] -> mat
  
  # complete values
  diag(mat)=0
  mat[lower.tri(mat)]  <- t(mat)[lower.tri(mat)]
  return(mat)
}

triangular=function(m, order) {
  m %>% 
    `row.names<-`(colnames(.)) %>% 
    .[ordered_labels,ordered_labels] -> m
  m[upper.tri(m)]=NA
  return(m)
}


#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------


#------------------------------------------------------------
# Function to convert list to dataframe
#------------------------------------------------------------
list2df <- function(L) {
  suppressPackageStartupMessages(require("reshape2", quietly=TRUE))
  param <- c("depth","lat","lon")
  n <- length(param)
  variable <- names(L)
  
  if (!is.list(L)) {
    stop("Need a list")
  }
  
  if (!all(param %in% variable)) {
    stop("Need a list with 4 components named lon, lat, depth and w, the variable of interest")
  }
  
  df <- melt(L[variable[4]], varnames=param)[,-5]
  df$lon <- L$lon[df$lon]
  df$lat <- L$lat[df$lat]
  df$depth <- L$depth[df$depth]
  colnames(df) <- c("depth","lat","lon",variable[4])
  
  return(df)
}


formating_name = function(x) {
  paste0("TARA_", 
         sapply(x, function(y) 
           paste0(paste0(rep(0,3-nchar(y)), collapse=""), y) 
         )
  )
}