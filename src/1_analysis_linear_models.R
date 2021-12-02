# Statistical analysis of protist population
# metagenomics
#
# Linear modeling
#
# December, 2021
#----------------------------------------

# options
options(stringsAsFactors = F)

# import libraries
library("tidyverse")
library("leaps")

# cleaning
rm(list=ls())
graphics.off()

# Read files
gen <- read.table("data/out/median_pairwise_fst_meanloc.txt", header=T)
env <- read.table("data/out/env_long.txt", header = T)
geo <- read.table("data/out/geo-distances.txt", header = T)
mct <- read.table("data/out/mct_df.txt", header = T)
pcnx <- read.table("data/out/proba_cnx.txt", header = T) 


# Formating pairwise station names
pairwiseName= function(df) {
  df %>% 
    filter(from != to) %>% 
    mutate(pStation = paste(gsub("TARA_","TO.",from),gsub("TARA_","TO.",to),sep="_")) %>% 
    select(-from, -to)
}

env <- pairwiseName(env)
geo <- pairwiseName(geo)
mct <- pairwiseName(mct)
pcnx <- pairwiseName(pcnx)
rm(pairwiseName)

# Joining all table
df <- Reduce(function(...) merge(..., all.x=TRUE, by = "pStation"), list(gen,geo,env,mct,pcnx))
rm(env, gen, geo, mct, pcnx)

# Processing NA
# Remove rows with missing Fst values (not enough material to compute Fst -> true missing values)
df <- df[-which(is.na(df$Fst)),]

# Replace by 0 for other factors (no difference between stations -> zero distance)
df[is.na(df)] = 0

write.table(df, "data/out/linear_model_df.txt", quote = F, row.names = F)


# Linear models
form = as.formula(Fst/(1-Fst) ~ geo + env + months_max + mct_min + 1)
form_null = as.formula(Fst/(1-Fst) ~ 1)

#----------------------
## Bathycoccus prasinos
#----------------------
df[which(df[,"Organism"] == "Bathycoccus-prasinos"),] %>%
  select(., -pStation, -Organism) %>%
  na.omit -> tmp

# Null versus complete model
complet.lm = lm(form, data = tmp)
null.lm = lm(form_null, data =  tmp)
anova(null.lm, complet.lm)

# Best model
tmp %>% 
regsubsets(form, ., method="exhaustive") %>%
  summary -> regfull
regfull$bic %>% which.min # 1
regfull # env

bathy.lm = lm(Fst/(1-Fst) ~ env, data = tmp)
anova(null.lm, bathy.lm)
summary(bathy.lm)

# Residuals vs variable
resid = data.frame(res = bathy.lm$residuals, env = tmp$env)
plot(resid$env,resid$res, xlab = "Distance environnementale", ylab = "Résidus")

# Predicts
x = seq(min(tmp$env),max(tmp$env),length.out = 1000)
pred = data.frame(envp = x,
                  Fstp = predict(bathy.lm, newdata = list(env=x)))


ggplot() +
  geom_point(aes(x = env, y = Fst/(1-Fst)), tmp) +
  geom_line(aes(envp, Fstp), pred, colour = "#F8766D") +
  xlab("Environmental distance") + ylab("Genetic distance")+
  ggtitle("Bathycoccus prasinos") +
  theme_classic() +
  theme(plot.title = element_text(face="bold.italic"))

rm(null.lm, bathy.lm, complet.lm)

#------------------------
## Pelagomonas calceolata
#------------------------
df[which(df[,"Organism"] == "Pelagomonas-calceolata"),] %>%
  select(., -pStation, -Organism) %>%
  na.omit -> tmp

complet.lm = lm(form, data = tmp)
null.lm = lm(form_null, data =  tmp)
anova(null.lm, complet.lm)

# Best model
tmp %>% 
  regsubsets(form, ., method="exhaustive") %>%
  summary -> regfull
regfull$bic %>% which.min # 1
regfull # geo

pelago.lm = lm(Fst/(1-Fst) ~ geo, data = tmp)
anova(null.lm, pelago.lm)
summary(pelago.lm)

rm(null.lm, pelago.lm, complet.lm)

#---------------------
## Phaeocystis cordata
#---------------------
df[which(df[,"Organism"] == "Phaeocystis-cordata"),] %>%
  select(., -pStation, -Organism) %>%
  na.omit -> tmp

# Null versus complete model
complet.lm = lm(form, data = tmp)
null.lm = lm(form_null, data =  tmp)
anova(null.lm, complet.lm)

# Model selection
tmp %>% 
  regsubsets(form, ., method="exhaustive") %>%
  summary -> regfull

regfull$bic %>% which.min # 1
regfull # geo

# Phaeocystis model
phaeo.lm = lm(Fst/(1-Fst) ~ geo, data = tmp)
anova(null.lm, phaeo.lm)
summary(phaeo.lm)

# Residuals vs variable
resid = data.frame(res = phaeo.lm$residuals, geo = tmp$geo)
plot(resid$geo,resid$res, xlab = "Distance géographique", ylab = "Résidus")

# Predicts
phaeo.lm = lm(Fst/(1-Fst) ~ geo, data = tmp)
x = seq(min(tmp$geo),max(tmp$geo),length.out = 1000)
pred = data.frame(geop = x,
                  Fstp = predict(phaeo.lm, newdata = list(geo=x)))
ggplot() +
  geom_point(aes(x = geo, y = Fst), tmp) +
  geom_line(aes(geop, Fstp), pred, colour = "#53B400") +
  xlab("Geographical distance") + ylab("Genetic distance")+
  ggtitle("Phaeocystis cordata") +
  theme_classic() +
  theme(plot.title = element_text(face="bold.italic"))

rm(null.lm, complet.lm)
