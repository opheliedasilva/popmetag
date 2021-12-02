#
# Probabilities of connection preprocessins
#
#
#
# Ophélie Da Silva
# December, 2021
#-------------------------------------------------------

# options
options(stringsAsFactors = F)

# import libraries
library("tidyverse")

# cleaning
rm(list=ls())
graphics.off()

# Initialisation
input = "data/raw/proba_cnx"
output = "data/out"

# Files
files = paste(input, list.files(input, pattern=".csv"), sep="/") %>%
  `names<-`(c("months3","months6","months12"))

# Read files:
df.l = list()
for(f in files) {
  read.csv(f, row.names = 1, header = T) %>%
    `colnames<-`(c("from","to", names(files[which(f == files)]))) -> df.l[[f]]
} ; rm(files,f,input)

# each matris is assymetric. For each pairwise station, only one value is kept
# The maximum porbabilities of connection are kept.
long.l <- list()
for (pc_df in df.l) {
  pc_param <- colnames(pc_df)[3] # current metric (pld: 3, 6, 12 months)

  pc_df %>% 
    # dataframe containing the two values for each pair of station
    `colnames<-`(c("to", "from", "x2")) %>% 
    left_join(pc_df %>% 
                `colnames<-`(c("from", "to", "x1")), 
              .) %>% 
    # identify the maximal value
    rowwise() %>%
    mutate(x_max = max(x1, x2, na.rm = T)) %>% 
    select(from, to, x_max) %>% 
    ungroup %>% 
    `colnames<-`(c("from", "to", 
                   paste0(pc_param,"_max"))) -> long.l[[pc_param]]
}

# homogenize station labels
label_format <- function(x) {
  # Function to format station label
  # from SXSUR to TARA_XXX
  x %>% 
    gsub("S([0-9]+)SUR", "\\1", .) %>%
    sapply(., function(x) {paste0(paste(rep(0,3-nchar(x)),collapse = ""),x)}) %>%
    paste0("TARA_",.)
}

long.l %>% 
  reduce(full_join) %>% 
  mutate(from = label_format(from),
         to = label_format(to)) -> pc.df

# TARA_027 and TARA_030 are very closed. They are in same cell grid. 
# But whereas TARA_027 is not a metaG station, TARA_030 is one. 
pc.df$from[which(pc.df$from == "TARA_027")] = "TARA_030"
pc.df$to[which(pc.df$to == "TARA_027")] = "TARA_030"

# averaged probability of connection over a year
pc.df %>% 
  replace(is.na(.), 0) %>% 
  mutate(months_max = (months3_max+months6_max+months12_max)/3) %>% 
  select(from, to, months_max) -> pc.df

# Writing 
write.table(pc.df, paste(output, "proba_cnx.txt", sep="/"), quote = F)
