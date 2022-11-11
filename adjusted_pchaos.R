#Computing series level 'adjusted prob of chaos'
#TR's attempt to reproduce

source('functions/functions_for_bayes_prob_.R')

S1 <- read.csv("Dataset S1.csv",skip = 1, na.strings = "NaN")
colnames(S1)[1] <- "MainID"
colnames(S1)[ncol(S1)] <- "adjpchaos"

prior.prob = 0.5
S1$adjpchaos2 = NA #reproduce results with incorrect LE results
S1$adjpchaos_correct = NA #use correct LE results
S1$adjpchaos_nost = NA #use only chaos tests, no 'stoch' results

# rates from all simulated datasets
LE = c(sensitivity = 1-0.29, specificity = 1-0.04, outcome = NA)
RQA = c(sensitivity = 1-0.37, specificity = 1-0.13, outcome = NA)
PE = c(sensitivity = 1-0.26, specificity = 1-0.18, outcome = NA)
# rates from tihelka (I think these include the periodic series)
#outcome = 0 if 'stochastic', omit test if 'determinsitic' or NA
# NPE = c(sensitivity = 1-0.05, specificity = 1-0.30, outcome = NA)
# KDE = c(sensitivity = 1-0.30, specificity = 1-0.13, outcome = NA)
# rates from tihelka with bug fix
NPE = c(sensitivity = 1-0.32, specificity = 1-0.28, outcome = NA)
KDE = c(sensitivity = 1-0.27, specificity = 1-0.58, outcome = NA)
# rates from simulated test dataset (first 10 reps, excluding periodic)
# NPE = c(sensitivity = 1-0.40, specificity = 1-0.55, outcome = NA)
# KDE = c(sensitivity = 1-0.175, specificity = 1-0.085, outcome = NA)

test.info.mat = data.frame(LE = LE, RQA = RQA, PE = PE, NPE = NPE, KDE= KDE)

for(i in 1:nrow(S1)) {
  test.info.mat.i <- test.info.mat.c <- test.info.mat
  test.info.mat.i[3,1:3] = ifelse(S1[i,c("LE1d","RQA","PE")]=="chaotic",1,0)
  test.info.mat.c[3,1:3] = ifelse(S1[i,c("LE","RQA","PE")]=="chaotic",1,0)
  test.info.mat.i[3,4:5] <- test.info.mat.c[3,4:5] <- ifelse(S1[i,c("NPE.class","KDE.class")]=="stochastic",0,NA)
  if(any(is.na(test.info.mat.i[3,]))) {
    test.info.mat.i = test.info.mat.i[,-which(is.na(test.info.mat.i[3,]))]
    test.info.mat.c = test.info.mat.c[,-which(is.na(test.info.mat.c[3,]))]
  }
  test.info.mat.n = test.info.mat.c[,1:3]
  S1$adjpchaos2[i] = get.posterior.prob(test.info.mat.i, prior.prob)
  S1$adjpchaos_correct[i] = get.posterior.prob(test.info.mat.c, prior.prob)
  S1$adjpchaos_nost[i] = get.posterior.prob(test.info.mat.n, prior.prob)
}

length(which(S1$adjpchaos2>0.5))/nrow(S1)
length(which(S1$adjpchaos_correct>0.5))/nrow(S1)
length(which(S1$adjpchaos_nost>0.5))/nrow(S1)

length(which(S1$adjpchaos2>0.9))/nrow(S1)
length(which(S1$adjpchaos_correct>0.9))/nrow(S1)
length(which(S1$adjpchaos_nost>0.9))/nrow(S1)
