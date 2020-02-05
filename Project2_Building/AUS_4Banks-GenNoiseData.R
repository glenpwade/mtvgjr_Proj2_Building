##====================== Initialisation ========================##
rm(list=ls())
gc(T)

if (TRUE) {
  library(graphics)
  library(Rcpp)
  library(RcppEigen)
  #library(moments)
  #library(Matrix)
  
  #setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC - GOOGLE DRIVE
  setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's laptop - GOOGLE DRIVE
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v4.r")
  source(functionsFile)
  
  cppfunctionsFile <- file.path(functionsPath,"cFunctions_tvgjr.cpp")
  Rcpp::sourceCpp(cppfunctionsFile, verbose=TRUE)
}

# Set the working directory - this is where the files will be saved.
setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        #Anna's Laptop - GOOGLE DRIVE

##====================== Initialisation ========================##


# We want iid noise to use for simulating the probability distributions
# Our samples contain 6149 observations.  We will match this and add 1500 more
# to provide a 'discard' dataset to stabilise the Garch process. (7659 in total)
# To keep things simple, we will round up the number of simulation observations to 8000
# We will use 5000 iterations in our simulations.  So we need 40,000,000 observations per sample.
# Keeping our 4 samples independant will be critical later, when we come to testing multivariate correlations.
# The best way to ensure this is to generate the data in one stream from a single seed.  We can then split
# our data into 4 discreet groups:

set.seed(1999)
noiseData <- matrix(rnorm(160000000),ncol = 20000,nrow = 8000)

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


##  Using the Garch Parameters we found by estimating a 'stable-looking' subset of the sample: sample[1201:2400]

##==============================================================##
##   ANZ
##==============================================================##
z <- readRDS("noiseANZ.RDS")
w <- z  # w will hold the data: Noise + Garch

#Pars: 0.1937064 0.123507 0.7934905 NA 
par_o <- 0.1937064
par_a <- 0.123507
par_b <- 0.7934905

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
  for (t in 1502:8000) {
    ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 
    w[t,n] <- sqrt(ht)*z[t,n]
  }
  
} #End: Generate Data

# Truncate the matrix to remove the discard rows & Save:
w <- w[1501:8000,]
saveRDS(w,"RefData_withGarch_ANZ_1.RDS")


##==============================================================##
##   CBA
##==============================================================##
z <- readRDS("noiseCBA.RDS")
w <- z  # w will hold the data: Noise + Garch

#Pars: 0.1454384 0.1080854 0.8039447 NA 
par_o <- 0.1454384
par_a <- 0.1080854
par_b <- 0.8039447

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
  for (t in 1502:8000) {
    ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 
    w[t,n] <- sqrt(ht)*z[t,n]
  }
  
} #End: Generate Data

# Truncate the matrix to remove the discard rows & Save:
w <- w[1501:8000,]
saveRDS(w,"RefData_withGarch_CBA_1.RDS")


##==============================================================##
##   NAB
##==============================================================##
z <- readRDS("noiseNAB.RDS")
w <- z  # w will hold the data: Noise + Garch

#Pars: 0.1857787 0.1529931 0.7674475 NA 
par_o <- 0.1857787
par_a <- 0.1529931
par_b <- 0.7674475

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
  for (t in 1502:8000) {
    ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 
    w[t,n] <- sqrt(ht)*z[t,n]
  }
  
} #End: Generate Data

# Truncate the matrix to remove the discard rows & Save:
w <- w[1501:8000,]
saveRDS(w,"RefData_withGarch_NAB_1.RDS")


##==============================================================##
##   WBC
##==============================================================##
z <- readRDS("noiseWBC.RDS")
w <- z  # w will hold the data: Noise + Garch

#Pars: 0.1244713 0.1275722 0.8101917 NA 
par_o <- 0.1244713
par_a <- 0.1275722
par_b <- 0.8101917

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
  for (t in 1502:8000) {
    ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 
    w[t,n] <- sqrt(ht)*z[t,n]
  }
  
} #End: Generate Data

# Truncate the matrix to remove the discard rows & Save:
w <- w[1501:8000,]
saveRDS(w,"RefData_withGarch_WBC_1.RDS")


