# Assess the accuracy of methods used for characterizing chaos and 
# stochasticity of time series


# Load libraries & set directory
if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, dplyr, tidyr, tibble, purrr, rpatrec, Chaos01, rdist, spam,  
               nonlinearTseries, stats, rdist, fields, viridis, viridisLite, readr)
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
                                                       "nonlienar_periodic")),
    maps = NA)     # empty column for simulated time series


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

# Run analyses on simulated series
source('functions/predict_np_udf.R')
source('functions/npe_heuristic_udf.R')

sim <- sim %>% 
  rowwise %>% 
  filter(!all(is.infinite(maps))) %>% # remove rows with infinite values
  mutate(ZeroOneTest = list(testChaos01(as.numeric(maps)))) %>% # 0-1 test for chaos
  mutate(NPE = list(npe_heuristic(maps))) # nonlinear prediction skill

sim$NPEclass=ifelse(sim$NPE>0.65,"deterministic","stochastic")
table(sim$NPEclass, sim$classification, useNA = "i")
#FPR comes from stoch series classified as det (which were all sine map, which is actually periodic), 
# excluding NAs (also only in sine map)

# Save output file
saveRDS(sim, "simulations_stoch.rds") # save as R data file

sim_exp <- sim %>% unnest(cols = c(maps, ZeroOneTest, NPE))
write_excel_csv(sim_exp,"simulations_stoch.csv") # save as .csv file
