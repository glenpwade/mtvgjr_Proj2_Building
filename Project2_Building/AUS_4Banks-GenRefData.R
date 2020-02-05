##====================== Initialisation ========================##

rm(list=ls())
gc(T)

#setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC & Laptop - GOOGLE DRIVE
#setwd("E:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Glen's Laptop - GOOGLE DRIVE
setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Glen's home PC - GOOGLE DRIVE


##====================== Initialisation ========================##


# We want iid noise to use for simulating the probability distributions
# Our samples contain 7027 observations.  We will match this and add 1500 more
# to provide a 'discard' dataset to stabilise the Garch process. (8527 in total)
# To keep things simple, we will round up the number of simulation observations to 10,000
# We will use 5000 iterations in our simulations.  So we need 50,000,000 observations per sample.
# Keeping our 4 samples independant will be critical later, when we come to testing multivariate correlations.
# The best way to ensure this is to generate the data in one stream from a single seed.  We can then split
# our data into 4 discreet groups:

if(FALSE) {
  
  set.seed(1999)
  noiseData <- matrix(rnorm(200000000),nrow = 10000,ncol = 20000)
  
  noiseSubset <- noiseData[,1:5000]
  saveRDS(noiseSubset,"noiseANZ.RDS")
  
  noiseSubset <- noiseData[,5001:10000]
  saveRDS(noiseSubset,"noiseCBA.RDS")
  
  noiseSubset <- noiseData[,10001:15000]
  saveRDS(noiseSubset,"noiseNAB.RDS")
  
  noiseSubset <- noiseData[,15001:20000]
  saveRDS(noiseSubset,"noiseWBC.RDS")
  
  
  ##  Tidy up and release variables  ##
  rm(noiseSubset,noiseData)
  gc(T)

}

##  Using the Garch Parameters we found by removing the rough-estimate of 'g':
##  Note: Rough-estimate calculated by smoothing the returns^2, then dividing the sqrt of this series
##  Out of the return series.  g' = r/sqrt(smooth(r^2))

##==============================================================##
##   ANZ
##==============================================================##
z <- readRDS("noiseANZ.RDS")
w <- z  # w will hold the data: Noise + Garch

#Pars: 
par_o <- 0.01744033
par_a <- 0.06944598
par_b <- 0.91295993

for (n in 1:5000) {
  # Initialise the first data points: e(t-1), h(t-1) using the first 1500 observations as discard
  ht_1 <- ht <- 1
  et_1 <- et <- z[1,n]
  for (t in 2:1500) {
    ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1
    et_1 <- et <- sqrt(ht)*z[t,n]
  }
  
  # Now generate the actual Garch data:  
  ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1 
  w[1501,n] <- sqrt(ht)*z[1501,n]
  for (t in 1502:10000) {
    ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 
    w[t,n] <- sqrt(ht)*z[t,n]
  }
  
} #End: Generate Data

# Truncate the matrix to remove the discard rows & Save:
w <- w[1501:10000,]
saveRDS(w,"RefData_withGarch_ANZ_2018.RDS")


##==============================================================##
##   CBA
##==============================================================##
z <- readRDS("noiseCBA.RDS")
w <- z  # w will hold the data: Noise + Garch

#Pars: 
par_o <- 0.02480881
par_a <- 0.08732556
par_b <- 0.88812449

for (n in 1:5000) {
  # Initialise the first data points: e(t-1), h(t-1) using the first 1500 observations as discard
  ht_1 <- ht <- 1
  et_1 <- et <- z[1,n]
  for (t in 2:1500) {
    ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1
    et_1 <- et <- sqrt(ht)*z[t,n]
  }
  
  # Now generate the actual Garch data:  
  ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1 
  w[1501,n] <- sqrt(ht)*z[1501,n]
  for (t in 1502:10000) {
    ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 
    w[t,n] <- sqrt(ht)*z[t,n]
  }
  
} #End: Generate Data

# Truncate the matrix to remove the discard rows & Save:
w <- w[1501:10000,]
saveRDS(w,"RefData_withGarch_CBA_2018.RDS")


##==============================================================##
##   NAB
##==============================================================##
z <- readRDS("noiseNAB.RDS")
w <- z  # w will hold the data: Noise + Garch

par_o <- 0.04194007
par_a <- 0.09714759
par_b <- 0.85929598

for (n in 1:5000) {
  # Initialise the first data points: e(t-1), h(t-1) using the first 1500 observations as discard
  ht_1 <- ht <- 1
  et_1 <- et <- z[1,n]
  for (t in 2:1500) {
    ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1
    et_1 <- et <- sqrt(ht)*z[t,n]
  }
  
  # Now generate the actual Garch data:  
  ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1 
  w[1501,n] <- sqrt(ht)*z[1501,n]
  for (t in 1502:10000) {
    ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 
    w[t,n] <- sqrt(ht)*z[t,n]
  }
  
} #End: Generate Data

# Truncate the matrix to remove the discard rows & Save:
w <- w[1501:10000,]
saveRDS(w,"RefData_withGarch_NAB_2018.RDS")


##==============================================================##
##   WBC
##==============================================================##
z <- readRDS("noiseWBC.RDS")
w <- z  # w will hold the data: Noise + Garch

#Pars: 
par_o <- 0.01931151
par_a <- 0.06944998
par_b <- 0.91135018

for (n in 1:5000) {
  # Initialise the first data points: e(t-1), h(t-1) using the first 1500 observations as discard
  ht_1 <- ht <- 1
  et_1 <- et <- z[1,n]
  for (t in 2:1500) {
    ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1
    et_1 <- et <- sqrt(ht)*z[t,n]
  }
  
  # Now generate the actual Garch data:  
  ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1 
  w[1501,n] <- sqrt(ht)*z[1501,n]
  for (t in 1502:10000) {
    ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 
    w[t,n] <- sqrt(ht)*z[t,n]
  }
  
} #End: Generate Data

# Truncate the matrix to remove the discard rows & Save:
w <- w[1501:10000,]
saveRDS(w,"RefData_withGarch_WBC_2018.RDS")


