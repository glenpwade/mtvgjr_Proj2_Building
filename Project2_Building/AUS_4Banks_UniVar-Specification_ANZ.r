rm(list=ls())

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


##====================== Initialisation ========================####
if (TRUE){
  library(graphics)
  library(foreach)
  library(doParallel)
  library(moments)
  library(Matrix)
  library(Rcpp)
  library(RcppEigen)
  
  setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC & Laptop - GOOGLE DRIVE

  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v3.r")
  source(functionsFile)
  
  cppfunctionsFile <- file.path(functionsPath,"cFunctions_tvgjr.cpp")
  Rcpp::sourceCpp(cppfunctionsFile, verbose=FALSE)
}
##====================== Initialisation ========================##

##====================== Data Setup ============================####
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
  mydata$e_cba <- mydata$e_cba - mean(mydata$e_cba)
  mydata$e_nab <- mydata$e_nab - mean(mydata$e_nab)
  mydata$e_wbc <- mydata$e_wbc - mean(mydata$e_wbc)
  
  # save data  
  saveRDS(mydata,"AUS_4Banks-ReturnsData.RDS")
  saveRDS(dates,"AUS_4Banks-Dates.RDS")
}
##====================== Data Setup ============================##

##====================== Data Load =============================####
if (TRUE){
  # Read AUS_4Banks data from saved file
  dates <- readRDS("AUS_4Banks-Dates.RDS")
  mydata <- readRDS("AUS_4Banks-ReturnsData.RDS")
  e_anz <- mydata$e_anz
  #e_cba <- mydata$e_cba
  #e_nab <- mydata$e_nab
  #e_wbc <- mydata$e_wbc
}
##====================== Data Load =============================##

##====================== Plot the Data =========================####
DoThis <- FALSE
if (DoThis) {
  ptitle <- "ANZ de-meaned Returns"
  plot(dates,e_anz,"l",main = ptitle)
  plot(e_anz,type="l",main = ptitle)
  ptitle <- "CBA de-meaned Returns"
  plot(dates,e_cba,"l",main = ptitle)
  ptitle <- "NAB de-meaned Returns"
  plot(dates,e_nab,"l",main = ptitle)
  ptitle <- "WBC de-meaned Returns"
  plot(dates,e_wbc,"l",main = ptitle)
}
##====================== Plot the Data =========================##



##==========================================================================##
##  This GenProbDist() function needs to be compiled before the code that follows it.  
##  It is just a subset of commands that needs to be run multiple times.   
##==========================================================================##

GenProbDist <- function(e,TV,refdata,testorder="H0",saveas="probDist",numloops=1100,UseCpp=FALSE) {
  
  # Calculate the test statistics:
  testOrder <- testorder
  refTestValues <- list()
  if(UseCpp){
    # Using Rcpp:
    refTestValues$LMTR2 <- Test_TV_noGARCH_TR2(TV$delta0,TV$pars,TV$shape,TV$st,e,TV$speedoption,testOrder)
    refTestValues$LMRobust <- Test_TV_noGARCH_ROBUST(TV$delta0,TV$pars,TV$shape,TV$st,e,TV$speedoption,testOrder)
    } else {
    refTestValues$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,testOrder)
    refTestValues$LMRobust <- myTest.TV.noGARCH.robust(e,TV,testOrder)
    }
   
  # ----- Set the parameters to pass to CalcProbabilityDist()  ------ #
  numCores <- 4
  numLoops <- numloops
  Tobs <- length(e)
  useTest <- c("TR2","ROBUST")    # This will force both Tests to be run!
  #useTest <- c("TR2")
  saveAs=paste0("Output/",saveas,".RDS")
  
  # ----- Load the generated data with Garch  ------ #
  RefData_WithGarch <- refdata 
  RefData_WithGarch <- RefData_WithGarch[1:Tobs,1:numLoops]
  ## Add 'g' into the Refernce Data.  
  RefData_WithGarch <- RefData_WithGarch*sqrt(TV$condvars)
  
  ### Generate the Probability Distribution and save to a file:
  timestamp(prefix = "##-- ",suffix = " --##")
  probDist <- CalcProbabilityDist_c(TV,RefData_WithGarch,refTestValues,useTest,testOrder,saveAs,numCores,numLoops) 
  timestamp(prefix = "##-- ",suffix = " --##")
  
  ## Extract Test P_Values from Results:
  Test <- list()
  Test$p_TR2 <- sum(probDist[,"Pval_TR2"],na.rm = TRUE)/(numLoops-length(which(is.na(probDist[,"Pval_TR2"]))))
  Test$p_ROB <- sum(probDist[,"Pval_Robust"],na.rm = TRUE)/(numLoops-length(which(is.na(probDist[,"Pval_Robust"]))))
 
  # Return:
  Test
  
}

##==============================================================##
##   General Specification Method:
##
##   1. Plot the local_Var and observe the transitions
##   2. Calculate the *most likely* garch parameters using different window lengths & kurtosis
##   3. Identify the first "single" transition start & end:
##     4. Generate refData0 using TV$delta0 = var(e[start:end])
##     5. Test for a first transition in this subset of the time series 
##     6. Generate refData1 using TV$delta0 & TV$pars = TV estimates found above
##     7. Test for a second transition in this subset of the time series 
##   8. Repeat steps 4 - 7 for each subsequent transition
##   9. Complete the specification by adding the individual sub-sections together
##  10. Consider if any of the transitions could be parameterised as doubles or squares?
##  11. Re-estimate the final model and calculate & check the standard errors

##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\ANZ ---- ORDER 0 ---- obs 1:6149           ##
##==============================================================##

DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(var(e),"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- NULL
  TV$shape <- vector("numeric")
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  parScale <- 1
  nDeps <- 1e-3
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  TEST <- GenProbDist(e,TV,refData,"H0",ptitle)

  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 0 ---- obs 1:6149           ##
## Execution time: 13 mins on Anna's laptop
## Result: Pvals=[TR2:19%,Robust=31%] 
## Fail to Reject the null: H0=Model is sufficient
##  This may be caused by issues with the optimiser where there are many transitions
##  So, try a smaller range...
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 0 ---- obs 1:1000           ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:1000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(var(e),"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- NULL
  TV$shape <- vector("numeric")
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  parScale <- 1
  nDeps <- 1e-3
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  TEST1_1000 <- GenProbDist(e,TV,refData,"H0",ptitle)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 0 ---- obs 1:1000            ##
## Result: Pvals=[TR2:26%, Robust:5.64%] 
## Fail to Reject the null: H0=Model is sufficient
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 0 ---- obs 1:3000           ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:3000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(var(e),"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- NULL
  TV$shape <- vector("numeric")
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  parScale <- 1
  nDeps <- 1e-3
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  TEST1_3000 <- GenProbDist(e,TV,refData,"H0",ptitle)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 0 ---- obs 1:3000            ##
## Result: Pvals=[TR2:32%, Robust:5.45%] 
## Fail to Reject the null: H0=Model is sufficient
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 0 ---- obs 1:4000           ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:4000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(var(e),"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- NULL
  TV$shape <- vector("numeric")
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  parScale <- 1
  nDeps <- 1e-3
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  TEST1_4000 <- GenProbDist(e,TV,refData,"H0",ptitle)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 0 ---- obs 1:4000            ##
## Result: Pvals=[TR2: 0.455%, Robust= 0%] 
## We strongly Reject the null: H0=Model is NOT sufficient
##
## Now we will test the alternate Hypothesis of:
## single, double & squared transition and see which offers the
## strongest level of rejection.
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 0, H03 ---- obs 1:4000           ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:4000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(var(e),"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- NULL
  TV$shape <- vector("numeric")
  TV$speedoption <- as.integer(3)
  TV$linpars <- c(NA,1,1,NA)
  TV$linparsMode <- 2
  parScale <- 1
  nDeps <- 1e-3
  TV$optimcontrol <- list(fnscale = -1) #, parscale=parScale, ndeps = nDeps)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H03 Test: Null H0=Double Transition, H1=Squared transition
  TEST1_4000_H03 <- GenProbDist(e,TV,refData,"H03",ptitle)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 0, H03 ---- obs 1:4000       
## Result: Pvals=[TR2:2.64%, Robust:4.36%] 
## Reject the null at 5% (but not 1%): Some Evidence of a triple
## transition, but this could also be a single...
##==============================================================##



##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 0, H02 ---- obs 1:4000           ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:4000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(var(e),"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- NULL
  TV$shape <- vector("numeric")
  TV$speedoption <- as.integer(3)
  TV$linpars <- c(NA,1,NA,NA)
  TV$linparsMode <- 2
  parScale <- 1
  nDeps <- 1e-4
  TV$optimcontrol <- list(fnscale = -1)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  if(is.null(refData)) refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H02 Test: Null H0=Single Transition, H1=Double transition
  TEST1_4000_H02 <- GenProbDist(e,TV,refData,"H02",ptitle)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 0, H02 ---- obs 1:4000       
## Result: Pvals=[TR2:68%, Robust:71%] 
## Fail to reject the null: No Evidence of a double transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 0, H01 ---- obs 1:4000           ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:4000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(var(e),"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- NULL
  TV$shape <- vector("numeric")
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  parScale <- 1
  nDeps <- 1e-3
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H01 Test: Null H0=No Transition, H1=single transition
  TEST1_4000_H01 <- GenProbDist(e,TV,refData,"H01",ptitle)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 0, H01 ---- obs 1:4000       
## Result: Pvals=[TR2:0%, Robust:0%] 
## Strongest Rejection of the null: Evidence for a single transition
##
## Now we will estimate the parameters for this transition, and
## then extend the sample window and repeat the process...
##==============================================================##

##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1 Estimation ---- obs 1:4000  ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:4000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.85  # Starting value selected from H01 Test run
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- c(-1,3,0.5)  #Starting values selected by observation of returns plot
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  parScale <- c(1,1,1,1)
  nDeps <- rep(1e-3,4)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  # Make a copy of the TV object.
  TV_1 <- TV
  
}  # End: if(DoThis)...
##==============================================================##
##   Estimated params are:
##   2.270099, -1.4617450, 6.9993616, 0.7158211
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1 ---- obs 1:4000           ##
##  Check for another transition in this sample subset:
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:4000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(2.270099,"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- c(-1.4617450,6.9993616,0.7158211,NA)
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
      
  if(is.null(refData)) refData <- readRDS("RefData_withGarch_ANZ.RDS")
  TEST1_4000_ORD1 <- GenProbDist(e,TV,refData,"H0",ptitle)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1 ---- obs 1:4000            ##
## Result: Pvals=[TR2: 84%, Robust= 46%] 
## Cannot Reject the null: H0=Model is sufficient
##
## Now we will extend the sample window and search for another
## transition.
##==============================================================##

##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1 ---- obs 1:5000           ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:5000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- c(2.0,3,0.5,NA)
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  parScale <- c(1,1,1,1)
  nDeps <- rep(1e-4,4)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  TV$optimcontrol <- list(fnscale = -1)

  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  TV <- TV_1
  
  if(!exists("refData")) refData <- readRDS("RefData_withGarch_ANZ.RDS")
  TEST5000_ORD1 <- GenProbDist(e,TV,refData,"H0",ptitle,2000)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1 ---- obs 1:5000            ##
## TestStat_ProbDist_ANZ Logliklihood Value:  -9169.499 
## Pars: 1.853966 3.866466 6.999961 0.8119224 NA 

## Result: Pvals=[TR2: 1.64%, Robust= 0%] 
## We Reject the null: H0=Model is NOT sufficient
## "Winner,Winner - Chicken Dinner!" Now, what shape is this fish?
## (25mins / 13.5mins using Rcpp version)
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1, H03 ---- obs 1:5000      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:5000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 2.0
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- c(0,5,0.5,NA)
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- c(NA,1,1,NA)
  TV$linparsMode <- 2
  parScale <- c(1,1,1,5,1,1)
  nDeps <- rep(1e-4,length(parScale))
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  TV$optimcontrol <- list(fnscale = -1) 
 
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars,TV$linpars, "\n\n")
  if (isTRUE(TV$error)) stop("There's been some mistake! - Estimation Failed!")
  
  TV <- TV_1

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H03 Test: Null H0=Double Transition, H1=Squared transition
  TEST5000_ORD1_H03 <- GenProbDist(e,TV,refData,"H03",ptitle,1100)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1, H03 ---- obs 1:5000       
## Result: Pvals=[TR2: 100%, Robust: 93%] 
#
## Cannot Reject the null - No Evidence of a triple transition
## 
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1, H02 ---- obs 1:5000      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:5000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(2.0,"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- c(-1,5,0.5,NA)
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- c(NA,1,NA,NA)
  TV$linparsMode <- 2
  #parScale <- c(1,1,1,1,1)
  parScale <- c(1,2,1,5,2)
  nDeps <- c(1e-5,1e-5,1e-3,1e-3,1e-5)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  TV$optimcontrol <- list(fnscale = -1)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars,TV$linpars, "\n\n")
  TV <- TV_1
  
  if(!exists("refData")) refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H02 Test: Null H0=Single Transition, H1=Double transition
  TEST5000_ORD1_H02 <- GenProbDist(e,TV,refData,"H02",ptitle,1200)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1, H02 ---- obs 1:5000       
## Result: Pvals=[TR2: 99.9%, Robust: 3.17%] 
#
## CANNOT Reject the null - No Evidence of a double
## transition, but the Robust test rejects at 5%
##==============================================================##

##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1, H01 ---- obs 1:5000      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:5000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- as.vector(2.0,"numeric")
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- c(1.0,3,0.7,NA)
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  parScale <- c(1,1,1,5)
  nDeps <- rep(1e-5,4)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  TV$optimcontrol <- list(fnscale = -1)
  
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  TV <- TV_1
  
  if(!exists("refData")) refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H01 Test: Null H0=No Transition, H1=Single transition
  TEST5000_ORD1_H01 <- GenProbDist(e,TV,refData,"H01",ptitle,1200)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1, H01 ---- obs 1:5000       
## Result: Pvals=[TR2: 0.33%, Robust: 0.0%] 
#
## We Reject the null - Strong Evidence of a single transition!
##==============================================================##

##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 2 Estimation ---- obs 1:5000  ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz[1:5000]
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 2.0
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- c(1,3,0.5,NA,1,3,0.75,NA) 
  TV$shape <- c(1,1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  parScale <- c(1,1,1,5,1,1,5)
  nDeps <- rep(1e-6,7)
  nDeps <- c(1e-5,1e-5,1e-4,1e-7,1e-5,1e-4,1e-7)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  TV$optimcontrol <- list(fnscale = -1)
  # Make a copy of the TV object.
  TV_1 <- TV
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  TV <- TV_1
  
}  # End: if(DoThis)...
##==============================================================##
##   Estimated params are:
##  Logliklihood Value:  -9169.501
##  Pars: 1.854318 0.07356302 6.976222 0.8119341 NA 3.784123 6.999993 0.8119445 NA 
##  NOTE:
##  This result is effectively the same as the Order 1 result:
##  Logliklihood Value:  -9169.499 
##  Pars: 1.853966 3.866466 6.999961 0.8119224 NA
##  So...
##  Lets expand the sample window and look again
##==============================================================##

##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1 ---- obs 1:6149           ##
##  Check for another transition in the full sample
##  using a TV Order 1 
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$st <- seq(1,Tobs)/Tobs
  TV$pars <- c(1,3,0.5,NA) 
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  parScale <- c(1,1,1,2)
  nDeps <- rep(1e-4,4)
  nDeps <- c(1e-5,1e-5,1e-4,1e-7)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  TV$optimcontrol <- list(fnscale = -1)
  # Make a copy of the TV object.
  TV_1 <- TV
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  TV <- TV_1

  if(!exists("refData")) refData <- readRDS("RefData_withGarch_ANZ.RDS")
  ptitle <- "TestStat_ProbDist_ANZ"
  TEST_ORD1 <- GenProbDist(e,TV,refData,"H0",ptitle,2000)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1 ---- obs 1:6149            ##
## Result: Pvals=[TR2: 0.25%, Robust: 0.0%] 
#
## We Reject the null - Strong Evidence of another transition!
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1, H03 ---- obs 1:6149      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.5,NA)
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- c(NA,1,1,NA)
  TV$linparsMode <- 2
  parScale <- c(1,1,1,1,1,1)
  nDeps <- rep(1e-4,length(parScale))
  TV$optimcontrol <- list(fnscale = -1, maxit=1500, parscale=parScale, ndeps=nDeps)
  TV$optimcontrol <- list(fnscale = -1, maxit=1500) 
 
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars,TV$linpars, "\n\n")
  TV <- TV_1

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H03 Test: Null H0=Double Transition, H1=Squared transition
  TEST_ORD1_H03 <- GenProbDist(e,TV,refData,"H03",ptitle,2000)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1, H03 ---- obs 1:6149       
## Result: Pvals=[TR2: 100%, Robust: 100%], with
## TestStat_ProbDist_ANZ Logliklihood Value:  -10995.17 
## Pars: 2.510853 4.155292 5.225108 0.6492135 NA NA 0.336503 -5.541565 NA 
#
## Cannot Reject the null - No Evidence of a triple transition
## 
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1, H02 ---- obs 1:6149      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.6,NA)
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- c(NA,1,NA,NA)
  TV$linparsMode <- 2
  parScale <- c(1,1,1,1,1)
  nDeps <- rep(1e-4,length(parScale))
  TV$optimcontrol <- list(fnscale = -1, maxit=1500, parscale=parScale, ndeps=nDeps)
  TV$optimcontrol <- list(fnscale = -1, maxit=1500) 
 
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars,TV$linpars, "\n\n")
  TV <- TV_1

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H02 Test: 
  TEST_ORD1_H02 <- GenProbDist(e,TV,refData,"H02",ptitle,2000)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1, H02 ---- obs 1:6149       
## Result: Pvals=[TR2: 100%, Robust: 100%], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -11036.6 
## Pars: 3.492153 3.268675 5.350097 0.6473324 NA NA -4.647861 NA NA 
#
## Cannot Reject the null - No Evidence of a double transition
## 
##==============================================================##

##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 1, H01 ---- obs 1:6149      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.5,NA)
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  parScale <- c(1,1,1,1)
  nDeps <- rep(1e-4,length(parScale))
  TV$optimcontrol <- list(fnscale = -1, maxit=1500, parscale=parScale, ndeps=nDeps)
  TV$optimcontrol <- list(fnscale = -1, maxit=1500) 
 
  # Make a copy of the TV object, before it is modified by the Estimation method
  TV_1 <- TV
  
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars,TV$linpars, "\n\n")
  TV <- TV_1

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H01 Test: 
  TEST_ORD1_H01 <- GenProbDist(e,TV,refData,"H01",ptitle,2000)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 1, H01 ---- obs 1:6149       
## Result: Pvals=[TR2: 0.80%, Robust: 0%], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -11259.21 
## Pars: 1.855147 1.53588 6.999802 0.6592161 NA
#
## Strongly Reject the null - Evidence of a single transition
## 
##==============================================================##


##==============================================================##
## TV - AUS_4Banks/ANZ ---- ORDER 2 Estimation ---- obs 1:6149  ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.6,NA,1,3,0.75,NA) 
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1,1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  TV$optimcontrol <- list(fnscale = -1)
  parScale <- c(1,1,1,5,1,1,5)
  nDeps <- rep(1e-3,7)
  nDeps <- c(1e-4,1e-4,1e-4,1e-3,1e-4,1e-4,1e-3)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  # Make a copy of the TV object.
  TV_1 <- TV
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  TV <- TV_1
  
}  # End: if(DoThis)...
##==============================================================##
##  ORDER 2 Estimated params are:
##  Logliklihood Value:  -11073.25
##  Pars: 2.271875 -1.506952 6.984017 0.465736 NA 2.637255 5.522067 0.647828 NA
##  NOW: Let's look for another transition...
##==============================================================##

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H0 Test: 
  TEST_ORD2 <- GenProbDist(e,TV,refData,"H0",ptitle,1500)

##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 2 ---- obs 1:6149       
## Result: Pvals=[TR2: 0.40%, Robust: 3.6%], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -11073.25 
## Pars: 2.271875 -1.506952 6.984017 0.465736 NA 2.637255 5.522067 0.647828 NA 
#
## Reject the null at 5% - Evidence of another transition
## So, now test for the shape...
##==============================================================##

##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 2, H03 ---- obs 1:6149      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.6,NA,1,3,0.75,NA) 
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1,1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- c(NA,1,1,NA)
  TV$linparsMode <- 2
  TV$optimcontrol <- list(fnscale = -1)

  parScale <- c(1,1,1,1,1,1,1,1,1)
  nDeps <- rep(1e-4,9)
  TV$optimcontrol <- list(fnscale = -1, maxit=2500, parscale=parScale, ndeps = nDeps)

  # Make a copy of the TV object.
  TV_1 <- TV
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars,TV$linpars, "\n\n")
  TV <- TV_1
  #plot(TV$condvars,type="l")

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H03 Test: 
  TEST_ORD2_H03 <- GenProbDist(e,TV,refData,"H03",ptitle,1500)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 2, H03 ---- obs 1:6149       
## Result: Pvals=[TR2: 100%, Robust: 85%], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -10891.56 
## Pars: 3.484879 -21.29173 2.619451 0.4555787 NA -53.26647 2.505412 0.8471711 NA NA -21.56073 89.88936 NA 
#
## Cannot Reject the null - No Evidence of a transition
## 
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 2, H02 ---- obs 1:6149      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.1,NA,1,3,0.7,NA) 
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1,1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- c(NA,1,NA,NA)
  TV$linparsMode <- 2
  TV$optimcontrol <- list(fnscale = -1)

  parScale <- c(1,1,1,1,1,1,1,1)
  nDeps <- rep(1e-4,8)
  #TV$optimcontrol <- list(fnscale = -1, maxit=2500, parscale=parScale, ndeps = nDeps)

  # Make a copy of the TV object.
  TV_1 <- TV
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars,TV$linpars, "\n\n")
  TV <- TV_1
  #plot(TV$condvars,type="l")

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H02 Test: 
  TEST_ORD2_H02 <- GenProbDist(e,TV,refData,"H02",ptitle,1500)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 2, H02 ---- obs 1:6149       
## Result: Pvals=[TR2: 100%, Robust: 100%], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -10891.56 
## Pars: 3.484879 -21.29173 2.619451 0.4555787 NA -53.26647 2.505412 0.8471711 NA NA -21.56073 89.88936 NA 
#
## Cannot Reject the null - No Evidence of a transition
## 
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks/ANZ ---- ORDER 2, H01 ---- obs 1:6149      ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.6,NA,1,3,0.75,NA) 
#  TV$pars <- c(1,3,0.1,NA,1,3,0.7,NA) 
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1,1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  TV$optimcontrol <- list(fnscale = -1)

  parScale <- c(1,1,1,5,1,1,5)
  nDeps <- c(1e-4,1e-4,1e-4,1e-3,1e-4,1e-4,1e-3)
  TV$optimcontrol <- list(fnscale = -1, maxit=2500, parscale=parScale, ndeps = nDeps)

  # Make a copy of the TV object.
  TV_1 <- TV
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars,TV$linpars, "\n\n")
  TV <- TV_1
  #plot(TV$condvars,type="l")

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H01 Test: 
  TEST_ORD2_H01 <- GenProbDist(e,TV,refData,"H01",ptitle,1500)
  
  ## Optional:  Look at the distributions
  probDist <- readRDS("Output/TestStat_ProbDist_ANZ.RDS")
  hist(probDist[,"Stat_TR2"])
  hist(probDist[,"Stat_Robust"])
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 2, H01 ---- obs 1:6149       
## Result: Pvals=[TR2: 14%, Robust: 2.9%], using
##  TestStat_ProbDist_ANZ Logliklihood Value:  -11073.25 
##  Pars: 2.271875 -1.506952 6.984017 0.465736 NA 2.637255 5.522067 0.647828 NA 
#
## Cannot Reject the null - No Evidence of a transition
## Can reject null at 5% using the Robust Test, so estimate Order 3:
##==============================================================##



##==============================================================##
## TV - AUS_4Banks/ANZ ---- ORDER 3 Estimation ---- obs 1:6149  ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.1,NA,1,3,0.5,NA,1,3,0.7,NA) 
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1,1,1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  #TV$optimcontrol <- list(fnscale = -1)

  parScale <- c(1,1,1,3,1,1,3,1,1,3)
  nDeps <- rep(1e-5,10)
  #nDeps <- c(1e-5,1e-5,1e-5,1e-6,1e-5,1e-5,1e-6,1e-5,1e-5,1e-6)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  # Make a copy of the TV object.
  TV_1 <- TV
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  TV <- TV_1

  plot(e,type="l")
  plot(TV$condvars,type="l")
  
}  # End: if(DoThis)...
##==============================================================##
##  ORDER 3 Estimated params are:
##  Logliklihood Value:  -10893.31
##  Pars: 1.900247 -0.04011741 5.914199 0.0001857999 NA 9.417206 6.996437 0.6613041 NA -9.627188 5.016682 0.7167887 NA 
##  NOW: Let's look for another transition...
##==============================================================##

  refData <- readRDS("RefData_withGarch_ANZ.RDS")
  # H0 Test: 
  TEST_ORD3 <- GenProbDist(e,TV,refData,"H0",ptitle,1500)

##==============================================================##
##  TV - AUS_4Banks/ANZ ---- ORDER 3 ---- obs 1:6149       
## Result: Pvals=[TR2: 7.2%, Robust: 0.0%], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -10893.31 
##  Pars: 1.900247 -0.04011741 5.914199 0.0001857999 NA 9.417206 6.996437 0.6613041 NA -9.627188 5.016682 0.7167887 NA 
#
## Cannot Reject the null at 5% with TR2, BUT
## Evidence of another transition given by Robust Test
## So, now test for the shape...
##==============================================================##


## Quick n dirty - just try to estimate Order 4 (single transitions)
##==============================================================##
## TV - AUS_4Banks/ANZ ---- ORDER 4 Estimation ---- obs 1:6149  ##
##==============================================================##
DoThis <- TRUE
if (DoThis) {
  
  e <- e_anz
  Tobs <- length(e)
  TV <- list()
  TV$delta0 <- 1.0
  TV$pars <- c(1,3,0.1,NA,1,3,0.3,NA,1,3,0.5,NA,1,3,0.7,NA) 
  TV$st <- seq(1,Tobs)/Tobs
  TV$shape <- c(1,1,1,1)
  TV$speedoption <- as.integer(3)
  TV$linpars <- NULL
  TV$linparsMode <- NULL
  TV$optimcontrol <- list(fnscale = -1)
  parScale <- rep(1,13)
  nDeps <- rep(1e-4,13)
  parScale <- c(1,1,1,3,1,1,3,1,1,3,1,1,3)
  #nDeps <- c(1e-5,1e-5,1e-5,1e-6,1e-5,1e-5,1e-6,1e-5,1e-5,1e-6)
  TV$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
  # Make a copy of the TV object.
  TV_1 <- TV
  ptitle <- "TestStat_ProbDist_ANZ"
  TV <- EstimateTV(e,TV,calcHess=FALSE)
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\nPars:",TV$delta0,TV$pars, "\n\n")
  TV <- TV_1

  plot(e,type="l")
  plot(TV$condvars,type="l")
  
}  # End: if(DoThis)...
## Result: Cannot get a higher likelihood value than with ORD 3

############  FINAL  SPECIFICATION  ############
##==============================================================##
##  ORDER 3 Estimated params are:
##  Logliklihood Value:  -10893.31
##  Pars: 1.900247 -0.04011741 5.914199 0.0001857999 NA 9.417206 6.996437 0.6613041 NA -9.627188 5.016682 0.7167887 NA 
##==============================================================##





##==============================================================##
##                            THE END
##==============================================================##




