# Testing GPDD time series for nonlinear stochasticity

# ***TR's comments

# Load libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rdist, tseriesChaos, fractal)
#The fractal package has been archived now and so you may have to download
#it manually from the CRAN archive

source('functions/embed_udf.R')
source('functions/predict_np_udf.R')

# List all time series
input_files <- list.files(path = "data/", pattern = "*.csv")

# ***store stuff
d <- m <- tw <- nse.x <- kap.prop <- numeric(length(input_files))

# Loop and run tests on each
for(i in seq_along(input_files)){ # ***run time ~3 mins
  df <- read.csv(paste0("data/", input_files[i]), na.strings = "N/A")
  
  #Step 1: Select time series to analyse
  x <- as.numeric(df$PopRescale)
  x <- na.omit(x) # ***note:omitting missing values leads to unequally spaced time steps
  
  print("**************************************")
  print(input_files[i])
  
  #Step 2: Embed time series (x)
  # ***add error catch if this fails
  results.embed_udf <- tryCatch(embed_udf(x),  error=function(err) NA) # delta-epsilon test
  if(all(is.na(results.embed_udf))) {
    d[i] <- m[i] <- tw[i] <- nse.x[i] <- NA
  } else {
    d[i]<-results.embed_udf[[1]] # ***time delay (tau)
    m[i]<-results.embed_udf[[2]] # ***embedding dimension (E)
    tw[i]<-results.embed_udf[[3]] #embedding dimension
    Mx<-results.embed_udf[[4]] #embedded data matrix
    dist<-rdist(Mx) #distance matrix, rdist(fields)
    
    #Step 3: Nonlinear prediction skill
    results.np<-predict_np_udf(Mx,m[i],frac.learn=0.5)
    nse.x[i]<-results.np[[1]]
  }
  print("nonlinear prediction skill:")
  print(nse.x[i]) #nonlinear prediction skill (Nash-Sutcliffe model efficiency)
  
  #Step 4: Kaplan's delta-epsilon test for determinism
  # ***add error catch if this fails, move to next i
  kap.test <- tryCatch(determinism(x),  error=function(err) NA) # delta-epsilon test
  if(all(is.na(kap.test))) {
    kap.prop[i] <- NA
    next
  }
  #print(kap.test) # summary of the analysis
  
  # ***try to get this without plotting
  kt <- plot(kap.test)
  Evals <- unlist(kt$E.orig.box) #observed
  Esurr <- unlist(lapply(kt$E.bxp, FUN = function(x) apply(x$stats,2,min))) #mins of surrogate data
  kap.prop[i] <- length(which(Evals<Esurr))/length(Evals) #prop observed less than min of surrogates across all conditions
  
  file_name = paste("result_plots/Kaplan_test_GPDD_", tools::file_path_sans_ext(input_files[i]), 
                    ".tiff", sep="") # name plots
  tiff(file_name, units="in", width=5, height=5, res=300, compression = "lzw")
  plot(kap.test) # plots E-statistics of  original series and surrogates
  dev.off()
}

# ***make output df
outdf <- data.frame(MainID = as.numeric(sub(".csv","",input_files)),
                    d = d, m = m, tw = tw, NPE = nse.x, kap.prop = kap.prop)
outdf <- outdf[order(outdf$MainID),]
outdf$NPEclass <- ifelse(outdf$NPE>0.65,"deterministic","stochastic")
outdf$KDEclass <- ifelse(outdf$kap.prop>0.3,"deterministic","stochastic")
write.csv(outdf,"reproduce_gpdd_results.csv", row.names = F)

# ***compare to existing results
S1 <- read.csv("Dataset S1.csv",skip = 1, na.strings = "NaN")
colnames(S1)[1] <- "MainID"
colnames(S1)[ncol(S1)] <- "adjpchaos"

outdf2 <- dplyr::left_join(outdf, S1, by="MainID")
outdf3 <- outdf2[,c("MainID","NPE.x","NPE.y","NPEclass","NPE.class","KDEclass","KDE.class","kap.prop")]
table(new=outdf2$NPEclass,original=outdf2$NPE.class,useNA = "i")
table(new=outdf2$KDEclass,original=outdf2$KDE.class,useNA = "i")
