# Validate nonlinear stochasticity methods on simulated data

# Load libraries
pacman::p_load(rdist, tseriesChaos, fractal, dplyr, tidyr, purrr, ggplot2)

source('functions/embed_udf.R')
source('functions/predict_np_udf.R')

# Simulated time series
sims=read.csv("./simdata/simulation_dataset_test.csv")

sims_d=select(sims, ID, Model, TimeStep, NoiseLevel, Classification, Sim.1:Sim.100) %>% 
  gather(SimNumber, Value, Sim.1:Sim.100) %>% 
  group_by(ID,Model,Classification,NoiseLevel,SimNumber) %>%  mutate(TSlength=length(Value)) %>% ungroup() %>% 
  group_by(ID,Model,Classification,NoiseLevel,TSlength,SimNumber) %>% nest() %>% 
  mutate(data=map(data, as.data.frame)) %>% ungroup()
modelorder=unique(arrange(sims_d, Classification, Model)$Model)
sims_d$Model=factor(sims_d$Model, levels=modelorder)

#use first 10 reps of each model (about 1 hr, non-parallel)
sims_test=filter(sims_d, SimNumber %in% paste0("Sim.",1:10))

# store stuff
d <- m <- tw <- nse.x <- kap.prop <- nrow(sims_test)

# Loop and run tests on each
for(i in 1:nrow(sims_test)){ 
  df <- sims_test$data[[i]]
  
  #Step 1: Select time series to analyse
  x <- as.numeric(df$Value)
  x <- na.omit(x) # ***note:omitting missing values leads to unequally spaced time steps
  
  print(i)
  
  #Step 2: Embed time series (x)
  # ***add error catch if this fails
  results.embed_udf <- tryCatch(embed_udf(x), error=function(err) NA) # delta-epsilon test
  if(all(is.na(results.embed_udf))) {
    d[i] <- m[i] <- tw[i] <- nse.x[i] <- NA
  } else {
    d[i]<-results.embed_udf[[1]] # ***time delay (tau)
    m[i]<-results.embed_udf[[2]] # ***embedding dimension (E)
    tw[i]<-results.embed_udf[[3]] #embedding dimension
    Mx<-results.embed_udf[[4]] #embedded data matrix
    dist<-rdist(Mx) #distance matrix, rdist(fields)
    
    #Step 3: Nonlinear prediction skill
    results.np<-tryCatch(predict_np_udf(Mx,m[i],frac.learn=0.5), error=function(err) NA)
    nse.x[i]<-results.np[[1]]
  }
  
  #Step 4: Kaplan's delta-epsilon test for determinism
  #add error catch if this fails, move to next i
  kap.test <- tryCatch(determinism(x),  error=function(err) NA) # delta-epsilon test
  if(all(is.na(kap.test))) {
    kap.prop[i] <- NA
    next
  }

  kt <- plot(kap.test)
  Evals <- unlist(kt$E.orig.box) #observed
  Esurr <- unlist(lapply(kt$E.bxp, FUN = function(x) apply(x$stats,2,min))) #mins of surrogate data
  kap.prop[i] <- length(which(Evals<Esurr))/length(Evals) #prop observed less than min of surrogates across all conditions
}

#make output df
outdf <- data.frame(d = d, m = m, tw = tw, NPE = nse.x, kap.prop = kap.prop)
outdf <- cbind(sims_test[,1:6], outdf)
outdf$NPEclass <- ifelse(outdf$NPE>0.65,"deterministic","stochastic")
outdf$KDEclass <- ifelse(outdf$kap.prop>0.3,"deterministic","stochastic")
outdf$Classification2=ifelse(outdf$Classification=="chaotic", "chaotic", "not chaotic")
outdf$NoiseLevel2=ifelse(outdf$NoiseLevel==0, 0.01, outdf$NoiseLevel)

write.csv(outdf,"sims_test_results.csv", row.names = F)

#overall prop correct classification ####
sims_long=gather(outdf, Method, Methodclass, NPEclass:KDEclass)
sims_long %>% 
  group_by(Method, Classification, Methodclass) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method, nesting(Classification, Methodclass), fill=list(n=0)) %>% 
  group_by(Method, Classification) %>% mutate(proportion=n/sum(n)) %>% as.data.frame()

#      Method Classification   Methodclass    n  proportion
# 1  KDEclass        chaotic deterministic 1155 0.825000000
# 2  KDEclass        chaotic    stochastic  245 0.175000000
# 3  KDEclass        chaotic          <NA>    0 0.000000000
# 4  KDEclass       periodic deterministic 1318 0.941428571
# 5  KDEclass       periodic    stochastic   82 0.058571429
# 6  KDEclass       periodic          <NA>    0 0.000000000
# 7  KDEclass     stochastic deterministic  119 0.085000000
# 8  KDEclass     stochastic    stochastic 1154 0.824285714
# 9  KDEclass     stochastic          <NA>  127 0.090714286
# 10 NPEclass        chaotic deterministic  823 0.587857143
# 11 NPEclass        chaotic    stochastic  563 0.402142857
# 12 NPEclass        chaotic          <NA>   14 0.010000000
# 13 NPEclass       periodic deterministic 1356 0.968571429
# 14 NPEclass       periodic    stochastic   10 0.007142857
# 15 NPEclass       periodic          <NA>   34 0.024285714
# 16 NPEclass     stochastic deterministic  767 0.547857143
# 17 NPEclass     stochastic    stochastic  593 0.423571429
# 18 NPEclass     stochastic          <NA>   40 0.028571429

#plot
summary=sims_long %>% 
  group_by(Method,Classification,Methodclass,NoiseLevel2,TSlength) %>% summarize(n=n()) %>% ungroup() %>% 
  complete(Method,nesting(Classification, Methodclass,NoiseLevel2,TSlength), fill=list(n=0)) %>% 
  group_by(NoiseLevel2,TSlength, Classification, Method) %>% 
  mutate(proportion=n/sum(n), Method2=sub("class","", Method))
ggplot(filter(summary, Methodclass=="stochastic"), aes(x=factor(NoiseLevel2), y=Classification, fill=proportion)) +
  facet_grid(Method2~TSlength) + geom_tile(stat = "identity") + 
  geom_text(aes(label=round(proportion,2)), color="white", size=3) +
  theme_classic() + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  labs(y="True Dynamics", x="Noise Level", fill="Proportion Classified Stochastic", title = "Time Series Length") +
  theme(plot.title = element_text(hjust = 0.5, size=11), strip.text = element_text(size=10), legend.position = "bottom",
        panel.spacing.x = unit(0.2,"lines"), panel.background = element_rect(color="black", fill=NULL),
        strip.background = element_blank()) 
