#
#
# Perform PCA and clustering
# Compute environmental distances
#
# Ophélie Da Silva
# December, 2021
# ---------------------------------------------------- 


# General seting
options(stringsAsFactors = F)
set.seed(0)

# Libraries
library("tidyverse")
library("bestNormalize")
library("FactoMineR")

# Cleaning
rm(list=ls())
graphics.off()

# Initial variables
env.f = "data/out/medar_env_tara_med.txt"
output = "data/out"

# read data
env <- read.table(env.f, header = T) %>% 
  select(-Depth)

# read gen data
"data/out/median_pairwise_fst_meanloc.txt" %>% 
  read.table(header=T)  %>% 
  mutate(from=str_split_fixed(pStation, "_", 2)[,1],
         to=str_split_fixed(pStation, "_", 2)[,2],
         from = gsub("TO\\.", "TARA_", from),
         to = gsub("TO\\.", "TARA_", to)) %>% 
  select(-pStation) -> fst

tara.station <- c(fst$from, fst$to) %>% unique
rm(fst)


# Format data
env %>% 
  # mutate(station = gsub("TARA_0+", "s", station)) %>% 
  select(station, longitude, latitude) -> tmp

# PCA - Environmental distances computing
# apply(env, 2, shapiro.test) # if p-value < .05, non normality
# # amon, cphl, ntra, ntri, slca, temp
# # psal, phos, dox1
#

var_to_norm <- c("amon", "cphl", "ntra", "ntri", "slca", "temp")

env %>%
  mutate_at(vars("psal","phos","dox1"), .funs = scale) %>%
  mutate_at(all_of(var_to_norm),
            .funs = function(x) {
              BN_obj <- bestNormalize(x, k = 3)
              return(predict(BN_obj))
              }
            ) -> env

# Perform PCA
env %>%
  `row.names<-`(.$station) %>% 
  PCA(., quali.sup = 1, quanti.sup = 2:3, graph = F) -> env.pca

# eigenvalues
pVar = env.pca$eig[, "percentage of variance"] ; m = mean(pVar)
sign = which(pVar > m)
sum(pVar[sign]) #82.67256

# coordinates of individuals along significant pc
env.pca$ind$coord[,sign] %>% 
  as.data.frame() %>%
  mutate(station = row.names(.)) %>%
  left_join(tmp, ., by = "station") %>%
  `colnames<-`(c("station","lon","lat","PC1","PC2","PC3")) -> coord.ind

env.pca$var$coord[,sign] %>% 
  as.data.frame() %>% 
  mutate(var = row.names(.)) %>% 
  rename(PC1 = Dim.1, PC2 = Dim.2, PC3 = Dim.3) -> coord.var

env.pca$eig %>% 
  .[sign,] %>% 
  as.data.frame() %>% 
  `names<-`(c("eigvalues","percentVar","cumulVal")) -> eig

rm(env, env.pca, m, sign, tmp)

# Clustering
coord.ind %>% 
  select(contains("PC")) %>%
  `row.names<-`(coord.ind$station) -> tmp

hc <- hclust(dist(tmp, diag = T), method = "ward.D2")

inertie <- sort(hc$height, decreasing = TRUE)
plot(inertie[1:20], type = "s", xlab = "Nombre de classes", ylab = "Inertie")
points(c(2, 3, 5, 6), inertie[c(2, 3, 5, 6)], 
       col = c("green3", "red3", "blue3"), cex = 2, lwd = 3)
bind_rows(cutree(hc, 2),
          cutree(hc, 3),
          cutree(hc, 5),
          cutree(hc, 6)) %>% 
  t() %>% 
  as.data.frame() %>% 
  `names<-`(c("clus2", "clus3", "clus5", "clus6")) %>% 
  mutate(station = row.names(.)) -> env_cluster

coord.ind %>% 
  full_join(env_cluster) %>% 
  mutate(geno = ifelse(station %in% tara.station, "Y", "N")) -> coord.ind


# write table
write.table(coord.ind, paste(output, "env_pca_ind.txt",sep="/"), quote = F)
write.table(coord.var, paste(output, "env_pca_var.txt",sep="/"), quote = F)
write.table(eig, paste(output, "env_pca_eig.txt",sep="/"), quote = F, row.names = F)

# Compute distances (weighted average)
expand.grid(coord.ind$station,coord.ind$station) %>%
  `colnames<-`(c("from","station")) %>%
  mutate_if(is.factor, as.character) -> env.long

env.long %>% 
  left_join(., coord.ind %>% 
              select(-lon, -lat)) %>% 
  rename(to = station, to_pc1 = PC1, to_pc2 = PC2, to_pc3 = PC3) %>% 
  left_join(., coord.ind %>% 
              select(-lon, -lat, -clus2, -clus3, -clus5, -clus6, -geno) %>% 
              rename(from = station, from_pc1 = PC1, from_pc2 = PC2, from_pc3 = PC3)) -> env.long


# compute distances
env.long %>% 
  mutate(denv = sqrt((from_pc1 - to_pc1)^2 +
                       (from_pc2 - to_pc2)^2 +
                       (from_pc3 - to_pc3)^2)) %>% 
  select("from","to","denv") -> env.long

rm(coord.ind, pVar)

# Remove duplicated values
env.long %>% 
  spread(key = to, value = denv) %>% 
  `row.names<-`(.$from) %>% 
  select(-from) -> env.m


env.m = env.m[, names(env.m) %>% 
                gsub("TARA_0+", "", .) %>% 
                as.numeric %>% 
                order
              ] # Order matrix columns
env.m = env.m[row.names(env.m) %>% 
                gsub("TARA_0+", "", .) %>% 
                as.numeric %>% 
                order,] # Order matrix rows

env.m[lower.tri(env.m)] = NA

env.df = env.m %>% 
  mutate(from = row.names(.)) %>% 
  gather(., key = "to", value = "env", -from) %>%
  na.omit

rm(env.long)

# Write data on disk
write.table(env.df, file = paste(output,"env_long.txt",sep="/"), row.names = F, quote = F)

