#
# Preprocessing the mean connection time data
#
#
#
# Ophélie Da Silva
# Decembre, 2021
#----------------------------------------------------

# options
options(stringsAsFactors = F)

# load libraries
library("tidyverse")

# cleaning
rm(list=ls())
graphics.off()

# Initialisation
input = "data/raw"
output = "data/out"

# input files
coord_file = "coord.txt"       # file containing geographic coordinates
mct_file = "mct_protist.txt"   # file containing mean connection time (from Berline et al., 2014)

# Reading files
coord_file %>% 
  paste(input, ., sep = "/") %>% 
  read.table(header = T, row.names = 1) -> coord.df

mct_file %>% 
  paste(input, ., sep = "/") %>% 
  read.table(header = F) %>%  
  `colnames<-`(row.names(coord.df)) %>%   # replace column names by the station name (same order)
  `row.names<-`(row.names(coord.df)) %>%  # replace row names by the station name (same order)
  `diag<-`(0) -> mct.m

rm(coord_file, mct_file)

# MCT matrix is assymetric as it corresponds to directed relationship
# Create a symmetric matrix corresponding to minimal pairwise mct value

# from wide to long format
mct.m %>%
  mutate(from = row.names(.)) %>% 
  gather(to, mct1, -from) -> tmp

#
tmp %>% 
  # dataframe containing the two values for each pair of station
  rename(mct2 = mct1, from = to, to = from) %>% 
  left_join(tmp, .) %>% 
  # identify the minimal value
  rowwise() %>%
  mutate(mct_min = min(mct1, mct2)) %>% 
  ungroup %>% 
  # keep columns of interest
  select(from, to, mct_min) -> mct.min

# Writing 
write.table(mct.min, paste(output, "mct_df.txt", sep = "/"), quote = F)
