rm(list=ls())
gc(T)

##==========================================================================================##
##  Note:  The Multivar Specification has not been rigorously tested at this stage, so...
##         This is a 'quick n dirty' estimation just to get indicative results.
##         We will complete the full specification & estimation later.  23-Feb-2018
##==========================================================================================##

##==========================================================================================##
##  Notes for use:
##
##    'DoThis <- TRUE/FALSE' is used to control the processing of distinct blocks of code 
##
##    The CalcProbabilityDist() can take some time to run. 
##    This is why the compiler & doParallel packages have been used.
##    Calls to this function may be commented out in the code below, 
##    since we only need to generate them once.  After that we read it back from a file.
##    To reuse the code with a new dataset, un-comment these lines.
##    
##==========================================================================================##


##====================== Initialisation ========================##
if (TRUE){
  library(graphics)
  library(foreach)
  library(doParallel)
  library(moments)
  library(Matrix)
  

  #setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC - GOOGLE DRIVE
  setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's laptop - GOOGLE DRIVE

  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v4_gw.r")
  source(functionsFile)
  
  library(Rcpp)
  library(RcppEigen)
  cppfunctionsFile <- file.path(functionsPath,"cFunctions_tvgjr.cpp")
  Rcpp::sourceCpp(cppfunctionsFile, verbose=TRUE)
}
##====================== Initialisation ========================##

##====================== Data Setup ============================##
DoThis <- FALSE
if (DoThis){
  mydata <- read.csv("Bank-Returns.csv",header=TRUE)
  dates <- as.Date(mydata$Date, format = "%d/%m/%Y")

  #Create scaled-up individtal data series vectors (Percentage Returns):
  mydata$e_anz <- mydata$anz * 100
  mydata$e_cba <- mydata$cba * 100
  mydata$e_nab <- mydata$nab * 100
  mydata$e_wbc <- mydata$wbc * 100
  
  #De-mean the returns data:
  mydata$e_anz <- mydata$e_anz - mean(mydata$e_anz)
  mydata$e_cba <- mydata$e_wbc - mean(mydata$e_cba)
  mydata$e_nab <- mydata$e_nab - mean(mydata$e_nab)
  mydata$e_wbc <- mydata$e_wbc - mean(mydata$e_wbc)
  
  # save data  
  saveRDS(mydata,"AUS_4Banks-ReturnsData.RDS")
  saveRDS(dates,"AUS_4Banks-Dates.RDS")
}
##====================== Data Setup ============================##

##====================== Data Load =============================##
if (TRUE){
  # Read AUS_4Banks data from saved file
  dates <- readRDS("AUS_4Banks-Dates.RDS")
  mydata <- readRDS("AUS_4Banks-ReturnsData.RDS")
  e_anz <- mydata$e_anz
  e_cba <- mydata$e_cba
  e_nab <- mydata$e_nab
  e_wbc <- mydata$e_wbc
}
##====================== Data Load =============================##

##====================== Plot the Data =========================##
DoThis <- FALSE
if (DoThis) {
  ptitle <- "WBC de-meaned Returns"
  plot(dates,e_anz,"l",main = ptitle)
  plot(e_anz,type="l",main = ptitle)
  ptitle <- "CBA de-meaned Returns"
  plot(dates,e_cba,"l",main = ptitle)
  ptitle <- "NAB de-meaned Returns"
  plot(dates,e_nab,"l",main = ptitle)
  ptitle <- "WBC de-meaned Returns"
  plot(dates,e_wbc,"l",main = ptitle)
  plot(e_wbc,type="l",main = ptitle)

}
##====================== Plot the Data =========================##



##==============================================================##
## Univariate TV Specifications are:

## ANZ:
##  ORDER 3 Estimated params are:
##  Logliklihood Value:  -10893.31
##  Pars: 1.900247 -0.04011741 5.914199 0.0001857999 NA 9.417206 6.996437 0.6613041 NA -9.627188 5.016682 0.7167887 NA
##  Shape <- c(1,1,1)

## CBA:
##  ORDER 2 Estimated params are:
##  Pars:  1.17046 + 8.834901*G(6.999959,0.6605106) - 8.807512*G(4.133976,0.7140517)
##  Shape <- (1,1)

## NAB:
##  ORDER 3 Estimated params are:
##  Logliklihood Value:   -10294.88 
##  Pars: 1.642661 -0.8646679 6.999913 0.4652876 NA 16.51015 4.759838 0.665668 NA -15.67746 4.292712 0.7048427 NA
##  Shape <- (1,1,1)

## WBC: 
##  ORDER 3 Estimated params are:
##  Logliklihood Value:  -10628.25 
##  Pars: 2.138931e-05 1.726756 6.381144 0.5416291 NA 5.388848 6.998857 0.6607741 NA -5.314586 4.253032 0.7236713 NA 
##  Shape <- c(3,1,1)  
##==============================================================##

##==============================================================##
## Step 1:  Create the 4 univariate TV objects
##==============================================================##

#  ANZ:
e <- e_anz
# Setup TV object with starting parameters
TV <- list()
TV$Tobs <- length(e)
TV$delta0 <- 1.0
vecpars <- c(1,3,0.1,1,3,0.5,1,3,0.7) 
TV$shape <- c(TVshape$linear,TVshape$linear,TVshape$linear)
TV$pars <- constructTVpars(vecpars,TV$shape)
TV$st <- seq(1,TV$Tobs)/TV$Tobs
TV$speedoption <- TVspeedopt$eta
TV$linpars <- NA
# Estimate the univar series:
#parScale <- c(1,1,1,3,1,1,3,1,1,3)
#nDeps <- rep(1e-5,length(parScale))
#TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps, reltol = 1e-8)
ptitle <- "TV_Univar_Estimate_ANZ"

TV$var_target <- TRUE

TV <- EstimateTV(e,TV,calcHess=F,verbose = T)
cat("\n", ptitle, "Logliklihood Value: ", TV$Estimated$value, "\nPars:",TV$Estimated$delta0,TV$Estimated$pars, "\n\n")
TV$Estimated$optimoutput
# Look at the plots:
plot(e,type="l")
plot(TV$condvars,type="l")

##  New GARCH Test  ##
e <- e_anz

GARCH <- list()
GARCH$type <- GARCHtype$constant
GARCH$pars <- var(e)

GARCH <- list()
GARCH$type <- GARCHtype$GJR_alpha0
GARCH$pars <- c(2.0,0.0,0.80,0.05)
#
# parScale <- c(1,1,5)
# nDeps <- c(1e-5,1e-5,1e-5)
# GARCH$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-8)
#
GARCH$var_target <- TRUE

g1 <- EstimateGARCH(e,GARCH,F,T)
a


#  CBA:
e <- e_cba
# Setup TV object with starting parameters
TV <- list()
TV$Tobs <- length(e)
TV$delta0 <- 1.0
TV$pars <- c(2,3,0.5,NA,-1,3,0.6,NA)
TV$shape <- c(1,1)
TV$st <- seq(1,TV$Tobs)/TV$Tobs
TV$speedoption <- 3
TV$linpars <- NULL
# Estimate the univar series:
parScale <- c(1,1,1,3,1,1,3)
nDeps <- rep(1e-4,length(parScale))
TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
ptitle <- "TV_Univar_Estimate_CBA"
TV <- EstimateTV(e,TV,calcHess=FALSE)
cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(TV$condvars,type="l")

TV_cba <- TV


#  NAB:
e <- e_nab
# Setup TV object with starting parameters
TV <- list()
TV$Tobs <- length(e)
TV$delta0 <- 1.0
TV$pars <- c(1,3,0.1,NA,1,3,0.5,NA,1,3,0.8,NA)  
TV$shape <- c(1,1,1)
TV$st <- seq(1,TV$Tobs)/TV$Tobs
TV$speedoption <- 3
TV$linpars <- NULL
# Estimate the univar series:
parScale <- c(1,1,1,1,1,1,1,1,1,1)
nDeps <- rep(1e-4,length(parScale))
TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
ptitle <- "TV_Univar_Estimate_NAB"
TV <- EstimateTV(e,TV,calcHess=FALSE)
cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(TV$condvars,type="l")

TV_nab <- TV


#  WBC:
e <- e_wbc
# Setup TV object with starting parameters
TV <- list()
TV$Tobs <- length(e)
TV$delta0 <- 1.0
TV$pars <- c(1,3,0.4,NA,1,3,0.5,NA,1,3,0.66,NA)
TV$shape <- c(3,1,1)
TV$st <- seq(1,TV$Tobs)/TV$Tobs
TV$speedoption <- 3
TV$linpars <- NULL
# Estimate the univar series:
parScale <- c(1,1,1,1,1,1,1,1,1,3)
nDeps <- rep(1e-4,length(parScale))
TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
ptitle <- "TV_Univar_Estimate_WBC"
TV <- EstimateTV(e,TV,calcHess=FALSE)
cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(TV$condvars,type="l")

TV_wbc <- TV


##==============================================================##
## Step 2:  Create the 4 univariate GARCH objects
##==============================================================##

# ANZ:

e <- e_anz/sqrt(TV_anz$condvars)
# Setup GARCH object with starting parameters
GARCH <- list()
GARCH$type <- 2  # GJR Garch
GARCH$pars <- c(0.05,0.05,0.85,0.05)  #Starting params
GARCH$var_target <- FALSE   # Omega (param #1) is free, not calculated based on alpha & beta
parScale <- c(3,3,1,3)
#parScale <- c(1,1,1,1)
nDeps <- rep(1e-4,length(parScale))
GARCH$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
# Estimate the univar series:
GARCH <- EstimateGARCH(e,GARCH,calcHess=FALSE)
ptitle <- "GARCH_Univar_Estimate_ANZ"
cat("\n", ptitle, "Logliklihood Value: ", GARCH$value, "\nPars:",GARCH$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(GARCH$condvars,type="l")

GARCH_anz <- GARCH


# CBA:

e <- e_cba/sqrt(TV_cba$condvars)
# Setup GARCH object with starting parameters
GARCH <- list()
GARCH$type <- 2  # GJR Garch
GARCH$pars <- c(0.05,0.05,0.85,0.05)  #Starting params
GARCH$var_target <- FALSE   # Omega (param #1) is free, not calculated based on alpha & beta
parScale <- c(3,3,1,3)
#parScale <- c(1,1,1,1)
nDeps <- rep(1e-4,length(parScale))
GARCH$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
# Estimate the univar series:
GARCH <- EstimateGARCH(e,GARCH,calcHess=FALSE)
ptitle <- "GARCH_Univar_Estimate_CBA"
cat("\n", ptitle, "Logliklihood Value: ", GARCH$value, "\nPars:",GARCH$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(GARCH$condvars,type="l")

GARCH_cba <- GARCH


# NAB:

e <- e_nab/sqrt(TV_nab$condvars)
# Setup GARCH object with starting parameters
GARCH <- list()
GARCH$type <- 2  # GJR Garch
GARCH$pars <- c(0.05,0.05,0.85,0.05)  #Starting params
GARCH$var_target <- FALSE   # Omega (param #1) is free, not calculated based on alpha & beta
parScale <- c(3,3,1,3)
#parScale <- c(1,1,1,1)
nDeps <- rep(1e-4,length(parScale))
GARCH$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
# Estimate the univar series:
GARCH <- EstimateGARCH(e,GARCH,calcHess=FALSE)
ptitle <- "GARCH_Univar_Estimate_NAB"
cat("\n", ptitle, "Logliklihood Value: ", GARCH$value, "\nPars:",GARCH$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(GARCH$condvars,type="l")

GARCH_nab <- GARCH


# WBC:

e <- e_wbc
# Setup GARCH object with starting parameters
GARCH <- list()
GARCH$type <- 2  # GJR Garch
GARCH$pars <- c(0.05,0.05,0.85,0.05)  #Starting params
GARCH$var_target <- FALSE   # Omega (param #1) is free, not calculated based on alpha & beta
parScale <- c(3,3,1,3)
#parScale <- c(1,1,1,1)
nDeps <- rep(1e-4,length(parScale))
GARCH$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
# Estimate the univar series:
GARCH <- EstimateGARCH(e,GARCH,calcHess=FALSE)
ptitle <- "GARCH_Univar_Estimate_WBC"
cat("\n", ptitle, "Logliklihood Value: ", GARCH$value, "\nPars:",GARCH$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(GARCH$condvars,type="l")

GARCH_wbc <- GARCH


##==============================================================##
## Step 3:  Create the multivariate objects for TV & GARCH
##==============================================================##


N <- 4  # We are modelling 4 Banks
nTV <- list(N)
nTV[[1]] <- TV_anz
nTV[[2]] <- TV_cba
nTV[[3]] <- TV_nab
nTV[[4]] <- TV_wbc
nTV$delta0 <- c(TV_anz$delta0,TV_cba$delta0,TV_nab$delta0,TV_wbc$delta0)


nGARCH <- list(N)
nGARCH[[1]] <- GARCH_anz
nGARCH[[2]] <- GARCH_cba
nGARCH[[3]] <- GARCH_nab
nGARCH[[4]] <- GARCH_wbc



##==============================================================##
## Step 4:  Create the STCC object
##==============================================================##

## Calculate the unconditional correlation matrix for standardized returns:
z_anz <- e_anz/(sqrt(TV_anz$condvars)*sqrt(GARCH_anz$condvars))
z_cba <- e_cba/(sqrt(TV_cba$condvars)*sqrt(GARCH_cba$condvars))
z_nab <- e_nab/(sqrt(TV_nab$condvars)*sqrt(GARCH_nab$condvars))
z_wbc <- e_wbc/(sqrt(TV_wbc$condvars)*sqrt(GARCH_wbc$condvars))

z <- cbind(z_anz,z_cba,z_nab,z_wbc)
uncondCorr <- cor(z)


STCC <- list()
STCC$Tobs <- length(e)                  # Note: Samples are constrained to be thesame length!
STCC$type <- 1                          # 0:"CCC", 1:"STCC",...
STCC$P1pars <- as.vector(uncondCorr)       # Identity Matrix - vectorised
STCC$P2pars <- as.vector(uncondCorr)       # Identity Matrix - vectorised
STCC$pars <- c(3,0.5,NA)                # A vector of the remaining pars: Speed, Loc1, <Loc2>
STCC$speedoption <- 3                   # speedoptions are: 1=gamma 2=gamma/sd(st) 3=exp(eta) 4=...
STCC$st <- seq(1,STCC$Tobs)/STCC$Tobs   # transition variable.  Default = seq(1,Tobs)/Tobs, i.e. linear
STCC$shape <- c(1)                      # vector of integers describing the transitions (G's)
#                                       # e.g. c(NA) => No transitions.  c(1,2,3) => 3 transitions
#                                       # 1=Single location, 2=Double location, 3=Single location, with Squared transition



##==============================================================##
## Step 5:  Estimate the STCC object
##==============================================================##






##==============================================================##
##                            THE END
##==============================================================##




