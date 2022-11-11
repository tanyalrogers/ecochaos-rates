# Assess the accuracy of methods used for characterizing chaos and 
# stochasticity of time series


# Load libraries & set directory
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, dplyr, tidyr, tibble, purrr, rpatrec, Chaos01, rdist, spam,  
               nonlinearTseries, stats, rdist, fields, viridis, viridisLite, readr, fractal)
# rpatrec version 1.0.1 doesn't work in some of the latest installments of R.
# This script was written and tested in R version 4.1.2 (Bird Hippie).


# Define variables
model <- c("logisticMap", "henonMap", "freitasMap", "sineMap", "cubicMap", 
           "poincareOscillator")
noise.level <- as.numeric(c("0", "0.001", "0.01", "1"))
length <- as.numeric(c("25", "50", "75", "100", "250"))


# Create data frame to hold simulations
sim <- tidyr::expand_grid(model, length, noise.level)
sim <- tibble::rowid_to_column(sim, "id")

sim <- sim %>%
  mutate(classification = model) %>%
  mutate(classification = case_when(
    classification %in% c("logisticMap", "henonMap")  ~ "chaotic",
    classification %in% c("freitasMap", "sineMap")  ~ "nonlinear_stochastic",
    classification %in% c("cubicMap", "poincareOscillator")  ~ "nonlienar_periodic",
    TRUE ~ NA_character_),   # classify map dynamics
    classification = factor(classification, levels = c("chaotic", 
                                                       "nonlinear_stochastic", 
                                                       "nonlienar_periodic")))


# Generate simulated series
source("functions/nonLinearMaps.R")

sim <- sim %>% 
  mutate(maps = invoke_map(model, length, noise.level)) %>% 
  unnest(maps) # note that some randomly generated initial values may lead to an
# unstable system that will tend to infinity

#plot the time series
par(mfrow=c(3,2))
test=filter(sim, length==25 & noise.level==0)
for(i in 1:nrow(test)) {
  plot(test$maps[[i]], type="o", ylab="value",
       main=paste(test$model[i],test$classification[i],"\nNoise = ",test$noise.level[i]))
}
test=filter(sim, length==25 & noise.level==1)
for(i in 1:nrow(test)) {
  plot(test$maps[[i]], type="o", ylab="value",
       main=paste(test$model[i],test$classification[i],"\nNoise = ",test$noise.level[i]))
}

#There is no observation noise being added. Also, there are transients, and the 
#sine map is actually a periodic 2-cycle with a transient, not stochastic.

#Actually adding the observation noise
sim2 <- sim %>% 
  rowwise %>% mutate(pars = list(list(n=length,noise.level=noise.level))) %>% ungroup() %>% 
  mutate(maps = invoke_map(model, pars),
         maps = map(maps, unlist)) %>% 
  rowwise %>% filter(!all(is.infinite(maps)))

#plot the time series
par(mfrow=c(3,2))
test=filter(sim2, length==25 & noise.level==0)
for(i in 1:nrow(test)) {
  plot(test$maps[[i]], type="o", ylab="value",
       main=paste(test$model[i],test$classification[i],"\nNoise = ",test$noise.level[i]))
}
test=filter(sim2, length==25 & noise.level==1)
for(i in 1:nrow(test)) {
  plot(test$maps[[i]], type="o", ylab="value",
       main=paste(test$model[i],test$classification[i],"\nNoise = ",test$noise.level[i]))
}

# Run analyses on simulated series
source('functions/predict_np_udf.R')
source('functions/npe_heuristic_udf.R')

kaplan=function(x) {
  kap.test <- tryCatch(determinism(x),  error=function(err) NA) # delta-epsilon test
  if(all(is.na(kap.test))) {
    kap.prop <- NA
  } else {
    kt <- plot(kap.test)
    Evals <- unlist(kt$E.orig.box) #observed
    Esurr <- unlist(lapply(kt$E.bxp, FUN = function(x) apply(x$stats,2,min))) #mins of surrogate data
    kap.prop <- length(which(Evals<Esurr))/length(Evals) #prop observed less than min of surrogates across all conditions
  }
  return(kap.prop)
}

sim <- sim %>% 
  rowwise %>% 
  filter(!all(is.infinite(maps))) %>% # remove rows with infinite values
  mutate(ZeroOneTest = list(testChaos01(as.numeric(maps)))) %>% # 0-1 test for chaos
  mutate(NPE = list(npe_heuristic(maps))) # nonlinear prediction skill
# This crashes R, possibly because there are many periodic series with no noise, 
# known to cause problem with off-the-shelf TDE.
# sim <- sim %>% 
#   rowwise %>% mutate(Kap.prop = kaplan(maps))

sim2 <- sim2 %>% 
  rowwise %>% 
  filter(!all(is.infinite(maps))) %>% # remove rows with infinite values
  mutate(ZeroOneTest = list(testChaos01(as.numeric(maps)))) %>% # 0-1 test for chaos
  mutate(NPE = list(npe_heuristic(maps))) %>%  # nonlinear prediction skill
  mutate(Kap.prop = kaplan(maps))

#without noise
sim$NPEclass=ifelse(sim$NPE>0.65,"deterministic","stochastic")
table(sim$NPEclass, sim$classification, useNA = "i")
#FPR comes from stoch series classified as det (which were all sine map, which is actually periodic), 
# excluding NAs (also only in sine map)
#FNR
length(which(sim$NPEclass=="stochastic" & sim$classification!="nonlinear_stochastic"))/length(which(sim$classification!="nonlinear_stochastic"))
#FPR
length(which(sim$NPEclass=="deterministic" & sim$classification=="nonlinear_stochastic"))/length(which(sim$classification=="nonlinear_stochastic"))
#numbers do not exactly match paper, likely because of random starting conditions with no seed set,
#insufficient number of replicates done for convergence.

#with the noise added
sim2$NPEclass=ifelse(sim2$NPE>0.65,"deterministic","stochastic")
sim2$KDEclass=ifelse(sim2$Kap.prop>0.3,"deterministic","stochastic")
table(sim2$NPEclass, sim2$classification, useNA = "i")
table(sim2$KDEclass, sim2$classification, useNA = "i")
#FNR
length(which(sim2$NPEclass=="stochastic" & sim2$classification!="nonlinear_stochastic"))/length(which(sim2$classification!="nonlinear_stochastic"))
#FPR
length(which(sim2$NPEclass=="deterministic" & sim2$classification=="nonlinear_stochastic"))/length(which(sim2$classification=="nonlinear_stochastic"))

#FNR
length(which(sim2$KDEclass=="stochastic" & sim2$classification!="nonlinear_stochastic"))/length(which(sim2$classification!="nonlinear_stochastic"))
#FPR
length(which(sim2$KDEclass=="deterministic" & sim2$classification=="nonlinear_stochastic"))/length(which(sim2$classification=="nonlinear_stochastic"))

# Save output file
saveRDS(sim, "simulations_stoch.rds") # save as R data file

sim_exp <- sim %>% unnest(cols = c(maps, ZeroOneTest, NPE))
write_excel_csv(sim_exp,"simulations_stoch.csv") # save as .csv file
