##====================== Initialisation ========================##

# We want to use a Wild Bootstrap for the NAB data

rm(list=ls())
gc(T)

if (TRUE){
  library(graphics)
  library(foreach)
  library(doParallel)
  library(moments)
  library(Matrix)
  
  #setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC - GOOGLE DRIVE
  #setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's laptop - GOOGLE DRIVE
  setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Glen's work PC - GOOGLE DRIVE
  
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v6.r")
  source(functionsFile)
  
    ## Setup the parallel backend environment ##
  numcores <- 3
  Sys.setenv("MC_CORES" = numcores)
  cl <- makeCluster(numcores)
  registerDoParallel(cl, cores = numcores)
  
  ## Remove parallel backend - if already setup ##
  if(FALSE) {
    unregisterDoParallel(cl)
    rm(cl)  
  }
  
  
}




##====================== Initialisation ========================##


#### ====================== Data Setup ============================ ####
DoThis <- FALSE
if (DoThis){
  mydata <- read.csv("Bank-Returns.csv",header=TRUE)
  dates <- as.Date(mydata$Date, format = "%d/%m/%Y")
  
  #Create scaled-up individtal data series vectors (Percentage Returns):
  mydata$e_nab <- mydata$anz * 100
  mydata$e_nab <- mydata$cba * 100
  mydata$e_nab <- mydata$nab * 100
  mydata$e_wbc <- mydata$wbc * 100
  
  #De-mean the returns data:
  mydata$e_anz <- mydata$e_anz - mean(mydata$e_anz)
  mydata$e_nab <- mydata$e_nab - mean(mydata$e_nab)
  mydata$e_nab <- mydata$e_nab - mean(mydata$e_nab)
  mydata$e_wbc <- mydata$e_wbc - mean(mydata$e_wbc)
  
  # save data  
  saveRDS(mydata,"AUS_4Banks-ReturnsData.RDS")
  saveRDS(dates,"AUS_4Banks-Dates.RDS")
}
##====================== Data Setup ============================##

#### ====================== Data Load ============================= ####
if (TRUE){
  # Read AUS_4Banks data from saved file
  dates <- readRDS("AUS_4Banks-Dates.RDS")
  mydata <- readRDS("AUS_4Banks-ReturnsData.RDS")
  e_nab <- mydata$e_nab
}
##====================== Data Load =============================##

#### ====================== Plot the Data ========================= ####
DoThis <- FALSE
if (DoThis) {
  ptitle <- "NAB de-meaned Returns"
  plot(dates,e_nab,"l",main = ptitle)
  ptitle <- "NAB de-meaned Returns"
  plot(e_nab,type="l",main = ptitle)
  abline(v=seq(1,6140,by = 100),lwd=0.25,col="lightgrey")
  lines(e_nab,type="l")
}
##====================== Plot the Data =========================##


##==============================================================##
##   General Specification Method:
##
##   1. Plot the demeaned returns and try to identify periods of constant volatility
##   2. Estimate the *most likely* sample garch parameters using a period from step 1
##   3. Iteratively test sub-sections of the sample for a transition (H0 Test)
##   4. Once a transition is found run the 3 aternate hypothesis tests to identify the *most likely* shape
##   5. Now re-estimate the sub-sample with the identified Order 1 model, and repeat the H0 Test
##   6  Continue refining the model until the sub-sample is fully specified
##   6. Repeat steps 4 - 6 for each subsequent sub-section
##   9. Complete the specification by adding the individual sub-sections together
##  10. Note: If the final specification is a high Order, it may be necessary to first
##       estimate it in 2 parts to identify *most likely* parameter estimates.
##       Then re-estimate the full model using tight-constraints in the optim controls.
##       <This is due to problems with optim in maximising high-order functions)
##  11. Consider if any of the transitions could/should be parameterised as squares or cubes?
##  12. Re-estimate the final model and finally calculate & check the standard errors

##==============================================================##


#### ========  SPECIFICATION - ORDER 0 1:3000  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 1:3000           ##
##==============================================================##

e <- e_nab[1:3000]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 1:3000            ##

## Result: Pvals=[TR2: 10.0%, Robust: 1.18%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -4958.437
## delta0: 1.596429
## Pars: NA
#
## Fail to Reject the null: H0=Model is sufficient, but...
##  The low Robust test score indicates evidence of a transition nearby
##  So, let's look at 1:2000 & 2000:4000...
##==============================================================##

##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 1:2000           ##
##==============================================================##
e <- e_nab[1:2000]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 1:2000            ##

## Result: Pvals=[TR2: 3.55%, Robust: 2.18%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -3184.451
## delta0: 1.413721
## Pars: NA
#
##  Reject the null: H0=Model is suspect
##  So, let's look for the shape...
##==============================================================## 



#### ---  H3,H2,H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H3 ---- obs 1:2000       ##
##==============================================================##

e <- e_nab[1:2000]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH03pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H03, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H3 ---- obs 1:2000         ##

## Result: Pvals=[TR2: 25.8%, Robust: 26%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -3130.408
## delta0: 0.9148737
## Pars: NA
## Linpars: -0.09732136  1.62535678
#
## Fail to Reject the NULL.  No evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H2 ---- obs 1:2000       ##
##==============================================================##

e <- e_nab[1:2000]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH02pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H2 ---- obs 1:2000         ##

## Result: Pvals=[TR2: 32.18%, Robust: 36.55%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -3133.77
## delta0: 0.7064782
## Pars: NA
## Linpars: 1.397693
#
## Cannot Reject the NULL.  No evidence for quadratic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H1 ---- obs 1:2000       ##
##==============================================================##

e <- e_nab[1:2000]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H01, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H1 ---- obs 1:2000         ##

## Result: Pvals=[TR2: 0.36%, Robust: 1.36%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -3184.451
## delta0: 1.413721
## Pars: NA
## Linpars: NA
#
## Reject the NULL.  Strongest evidence for linear transition
## Now, let's estimate this transition...
##==============================================================##


#### ========  SPECIFICATION - ORDER 1 1:2000  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 1:2000           ##
##==============================================================##

e <- e_nab[1:2000]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,5,0.7)
TV <- newTV(del0=1,vecpars=startpars,shape=TRshape$linear,speedopt=TRspeedopt$eta,st=ST)
# Estimate the univar series:
parScale <- c(1,1,1,1)
nDeps <- rep(1e-7,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-12)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l')

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 1:2000            ##

## Result: Pvals=[TR2: 80.58%, Robust: 62.88%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -3115.823
## delta0: 1.059839
## Pars: 
##           [,1]
## delta 1.2635942
## speed 6.9999999
## loc1  0.7220886
## loc2         NA
#
## Cannot Reject the null: H0=Model is sufficient
##==============================================================##


#### ========  SPECIFICATION - ORDER 0 2001:4000  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 2001:4000        ##
##==============================================================##

e <- e_nab[2001:4000]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 1100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 2001:4000         ##

## Result: Pvals=[TR2: 22.8%, Robust: 3.63%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -3204.643
## delta0: 1.443777
## Pars: NA
#
##  Mixed results from the Tests, so let's try a different sub-sample
##==============================================================## 

##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 2001:3000        ##
##==============================================================##

e <- e_nab[2001:3600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Null: Taylor order 3
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H0, saveas = fName,numloops = 5000)
beep()
# Null: Taylor order 1
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H01, saveas = fName,numloops = 1100)
beep()

# H01, 5000 loops:
# 2001:3500 = 6.32,0.74
# 2001:3600 = 5.26,0.8
# 2001:3700 = 7.04,1.9
# 2001:3800 = 6.28,2.18


##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 2001:3600         ##

## Result: Pvals=[TR2: 5.26%, Robust: 0.8%], using
## Taylor order 1 test.
## TestStat_ProbDist_NAB Logliklihood Value:  -2589.951
## delta0: 1.492101
## Pars: NA
#
##  Fail to reject the NULL: Evidence of a transition
##  Pretty sure it's linear, but let's test anyway...
##==============================================================## 

#### --- H3/H2/H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H3 ---- obs 2001:3600    ##
##==============================================================##

e <- e_nab[2001:3600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH03pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H03, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H3 ---- obs 2001:3600      ##

## Result: Pvals=[TR2: 42.91%, Robust: 36.27%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2483.05
## delta0: 3.415905
## Pars: NA
## Linpars: -5.682875  2.853870
#
## Fail to Reject the NULL.  No evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H2 ---- obs 2001:3600    ##
##==============================================================##

e <- e_nab[2001:3600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH02pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H2 ---- obs 2001:3600      ##

## Result: Pvals=[TR2: 48.2%, Robust: 32.3%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2491.858
## delta0: 2.48731
## Pars: NA
## Linpars: -2.058328
#
## Cannot Reject the NULL.  No evidence for quadratic transition
## Therefore it must be a Linear transition.
##==============================================================##


#### ========  SPECIFICATION - ORDER 1 2001:3600  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 2001:3600        ##
##==============================================================##

e <- e_nab[2001:3600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,4,0.5)
TV <- newTV(del0=1,vecpars=startpars,shape=TRshape$linear,speedopt=TRspeedopt$eta,st=ST)
# Estimate the univar series:
parScale <- c(1,1,1,3)
nDeps <- rep(1e-4,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-12)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l')

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Null Hyp: Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H0,saveas=fName,numloops = 1100)
# Null Hyp: Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H01,saveas=fName,numloops = 1100)


##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 2001:3600         ##

## Result: Pvals=[TR2: 99.8%, Robust: 99.7%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2468.056
## delta0: 2.183613
## Pars: 
##           [,1]
## delta -1.486306
## speed  6.999914
## loc1   0.538707
## loc2         NA
#
## Cannot Reject the null: H0=Model is sufficient
## So, let's look at the next sub-sample...
##==============================================================##



#### ========  SPECIFICATION - ORDER 0 3601:4600  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 3601:4600        ##
##==============================================================##

e <- e_nab[3601:4600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 1100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 3601:4600         ##

## Result: Pvals=[TR2: 1.0%, Robust: 0.0%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2243.59
## delta0: 5.205647
## Pars: NA
#
## Strongly Reject the Null: Evidence of a transition
## Let's check for the shape...
##==============================================================## 

#### --- H3/H2/H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H3 ---- obs 3601:4600    ##
##==============================================================##

e <- e_nab[3601:4600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH03pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H03, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H3 ---- obs 3601:4600      ##
##
## Result: Pvals=[TR2: 1.67%, Robust: 0.0%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2139.521
## delta0: 1.168338
## Pars: NA
## Linpars: -3.779833 21.522103
#
## Reject the NULL.  Evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H2 ---- obs 3601:4600    ##
##==============================================================##

e <- e_nab[3601:4600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH02pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H2 ---- obs 3601:4600      ##
##
## Result: Pvals=[TR2: 85.3%, Robust: 55.9%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2144.524
## delta0: 0.5588002
## Pars: NA
## Linpars: 9.176112
#
## Cannot Reject the NULL.  No evidence for quadratic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H1 ---- obs 3601:4600    ##
##==============================================================##

e <- e_nab[3601:4600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H01, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H1 ---- obs 3601:4600      ##
##
## Result: Pvals=[TR2: 5.18%, Robust: 0.0%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2243.59
## delta0: 5.205647
## Pars: NA
## Linpars: NA
##
## Reject the NULL.  Evidence for linear transition, but strongest
## support is for a cubic, so... estimation time...
## Or we could split this sub sample more & try to identify linear
## transitions.
##==============================================================##



#### ========  SPECIFICATION - ORDER 1 3601:4600  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 3601:4600        ##
##==============================================================##

e <- e_nab[3601:4600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,3,0.1)
TV <- newTV(del0=1,vecpars=startpars,shape=TRshape$cubic,speedopt=TRspeedopt$eta,st=ST)
# Estimate the univar series:
parScale <- c(1,1,1,3)
nDeps <- rep(1e-4,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-12)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l')

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Null Hyp: Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H0,saveas=fName,numloops = 1100)
beep()
# Null Hyp: Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H01,saveas=fName,numloops = 1100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 3601:4600         ##

## Result: Pvals=[TR2: 95.9%, Robust: 92.8%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2025.118
## delta0: 23.91409
## Pars: 
##           [,1]
## delta -22.8461455
## speed   3.5887420
## loc1    0.6570495
## loc2           NA
#
## Cannot Reject the null: H0=Model is sufficient
## So, let's look at the next sub-sample...
##==============================================================##


#### ========  SPECIFICATION - ORDER 0 4601:6149  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 4601:6149        ##
##==============================================================##

e <- e_nab[4601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 1100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 4601:6149         ##
##
## Result: Pvals=[TR2: 2.27%, Robust: 0.091%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2591.034
## delta0: 1.662163
## Pars: NA
#
## Strongly Reject the Null: Evidence of a transition
## Let's check for the shape...
##==============================================================## 


#### --- H3/H2/H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H3 ---- obs 4601:6149    ##
##==============================================================##

e <- e_nab[4601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH03pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H03, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H3 ---- obs 4601:6149      ##
##
## Result: Pvals=[TR2: 3.82%, Robust: 1.0%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2541.611
## delta0: 3.792863
## Pars: NA
## Linpars: -9.870255  8.523447
#
## Reject the NULL.  Evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H2 ---- obs 4601:6149    ##
##==============================================================##

e <- e_nab[4601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH02pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H2 ---- obs 4601:6149      ##
##
## Result: Pvals=[TR2: 0.18%, Robust: 0.27%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2589.913
## delta0: 1.782603
## Pars: NA
## Linpars: -0.2444253
#
## Reject the NULL.  Evidence for quadratic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H1 ---- obs 4601:6149    ##
##==============================================================##

e <- e_nab[4601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H01, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H1 ---- obs 4601:6149      ##
##
## Result: Pvals=[TR2: 55.8%, Robust: 62.8%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -2591.034
## delta0: 1.662163
## Pars: NA
## Linpars: NA
##
## Fail to Reject the NULL.  No Evidence for linear transition
## Strongest support is for a quadratic.
## so... estimation time...
##==============================================================##


#### ========  SPECIFICATION - ORDER 1 4601:6149  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 4601:6149        ##
##==============================================================##

e <- e_nab[4601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,3,0.1,0.5)
TV <- newTV(del0=1,vecpars=startpars,shape=TRshape$quadratic,speedopt=TRspeedopt$eta,st=ST)
# Estimate the univar series:
parScale <- c(1,1,1,3,3)
nDeps <- rep(1e-4,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-12)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l')

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Null Hyp: Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H0,saveas=fName,numloops = 1100)
beep()
# Null Hyp: Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H01,saveas=fName,numloops = 1100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 4601:6149         ##
##
## Tests fail to compute => Function is mis-specified.
## So, let's split this sub-sample...
##==============================================================##


#### ========  SPECIFICATION - ORDER 0 4601:5100  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 4601:5100        ##
##==============================================================##

e <- e_nab[4601:5100]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 2100)
beep()
# Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H01, saveas = fName,numloops = 2100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 4601:5100         ##
##
## Result: Pvals=[TR2: 4.38%, Robust: 1.19%], using Taylor Order 3
## TestStat_ProbDist_NAB Logliklihood Value:  -921.9145
## delta0: 2.342668
## Pars: NA
#
## Strongly Reject the Null: Evidence of a transition
## Let's check for the shape...
##==============================================================## 


#### --- H3/H2/H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H3 ---- obs 4601:5100    ##
##==============================================================##

e <- e_nab[4601:5100]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH03pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H03, saveas = fName,numloops = 2100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H3 ---- obs 4601:5100      ##
##
## Result: Pvals=[TR2: 64.6%, Robust: 59.9%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -905.2276
## delta0: 4.431782
## Pars: NA
## Linpars: -15.67047  17.90756
#
## Fail to Reject the NULL.  NO Evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H2 ---- obs 4601:5100    ##
##==============================================================##

e <- e_nab[4601:5100]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH02pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 2100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H2 ---- obs 4601:5100      ##
##
## Result: Pvals=[TR2: 24.3%, Robust: 12.7%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -919.9634
## delta0: 1.8983
## Pars: NA
## Linpars: 0.8708208
#
## Fail to Reject the NULL.  No Evidence for quadratic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H1 ---- obs 4601:5100    ##
##==============================================================##

e <- e_nab[4601:5100]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H01, saveas = fName,numloops = 2100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H1 ---- obs 4601:5100      ##
##
## Result: Pvals=[TR2: 43.4%, Robust: 48.4%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -921.9145
## delta0: 2.342668
## Pars: NA
## Linpars: NA
##
## Fail to Reject the NULL.  No Evidence for linear transition
## Strongest support is for a quadratic, though this failed to
## reject the null.
## so... estimation time...
##==============================================================##

#### ========  SPECIFICATION - ORDER 1 4601:5100  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 4601:5100        ##
##==============================================================##

e <- e_nab[4601:5100]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,3,0.3,0.6)
TV <- newTV(del0=1,vecpars=startpars,shape=TRshape$quadratic,speedopt=TRspeedopt$eta,st=ST)
# Estimate the univar series:
parScale <- c(1,1,1,1,1)
nDeps <- rep(1e-3,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-12)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l')

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Null Hyp: Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H0,saveas=fName,numloops = 1100)
beep()
# Null Hyp: Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H01,saveas=fName,numloops = 1100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 4601:5100         ##
##
## Tests fail to compute => Function is mis-specified.
## So, let's split this sub-sample...
##==============================================================##



#### ========  SPECIFICATION - ORDER 0 4601:5600  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 4601:5600        ##
##==============================================================##

e <- e_nab[4601:5600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 2100)
beep()
# Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H01, saveas = fName,numloops = 5000)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 4601:5600         ##
##
## Result: Pvals=[TR2: 27.7%, Robust: 14.7%], using Taylor Order 3
## Result: Pvals=[TR2: 6.1%, Robust: 1.74%], using Taylor Order 1
## TestStat_ProbDist_NAB Logliklihood Value:  -1682.798
## delta0: 1.696597
## Pars: NA
#
##  Borderline Reject the Null: Some Evidence of a transition
## Let's check for the shape...
##==============================================================## 


#### --- H3/H2/H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H3 ---- obs 4601:5600    ##
##==============================================================##

e <- e_nab[4601:5600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH03pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H03, saveas = fName,numloops = 2100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H3 ---- obs 4601:5600      ##
##
## Result: Pvals=[TR2: 47%, Robust: 51.8%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -1661.421
## delta0: 2.456319
## Pars: NA
## Linpars: -1.3326691 -0.2854366
#
## Fail to Reject the NULL.  NO Evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H2 ---- obs 4601:5600    ##
##==============================================================##

e <- e_nab[4601:5600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH02pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 2100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H2 ---- obs 4601:5600      ##
##
## Result: Pvals=[TR2: 93%, Robust: 92.8%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -1661.462
## delta0: 2.531706
## Pars: NA
## Linpars: -1.668446
#
## Fail to Reject the NULL.  No Evidence for quadratic transition
## So, if it exists, it must be linear...
##==============================================================##

#### ========  SPECIFICATION - ORDER 1 4601:5600  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 4601:5600        ##
##==============================================================##

e <- e_nab[4601:5600]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,3,0.5)
TV <- newTV(del0=1,vecpars=startpars,shape=TRshape$linear,speedopt=TRspeedopt$eta,st=ST)
# Estimate the univar series:
parScale <- c(1,1,1,3)
nDeps <- rep(1e-4,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-12)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l')

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Null Hyp: Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H0,saveas=fName,numloops = 1100)
beep()
# Null Hyp: Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H01,saveas=fName,numloops = 1100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 4601:5600         ##

## Result: Pvals=[TR2: 97.7%, Robust: 92.5%], using Taylor Order 3
## TestStat_ProbDist_NAB Logliklihood Value:  -1643.502
## delta0: 2.38904
## Pars: 
##           [,1]
## delta  -1.334988
## speed   4.534029
## loc1    0.480660
## loc2           NA
#
## Cannot Reject the null: H0=Model is sufficient
## So, let's look at the next sub-sample...
##==============================================================##



#### ========  SPECIFICATION - ORDER 0 5601:6149  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 5601:6149        ##
##==============================================================##

e <- e_nab[5601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,saveas = fName,numloops = 2100)
beep()
# Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H01, saveas = fName,numloops = 5000)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 5601:6149         ##
##
## Result: Pvals=[TR2: 1.95%, Robust: 0.10%], using Taylor Order 3
## TestStat_ProbDist_NAB Logliklihood Value:  -907.94
## delta0: 1.59947
## Pars: NA
#
##  Reject the Null: Evidence of a transition
## Let's check for the shape...
##==============================================================## 


#### --- H3/H2/H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H3 ---- obs 5601:6149    ##
##==============================================================##

e <- e_nab[5601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH03pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H03, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H3 ---- obs 5601:6149      ##
##
## Result: Pvals=[TR2: 91.8%, Robust: 92.4%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -854.3698
## delta0: 0.579344
## Pars: NA
## Linpars: -0.4521074  3.7313671
#
## Fail to Reject the NULL.  NO Evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H2 ---- obs 5601:6149    ##
##==============================================================##

e <- e_nab[5601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- addH02pars(TV)
TV <- EstimateTV(e,TV)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$linpars

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H2 ---- obs 5601:6149      ##
##
## Result: Pvals=[TR2: 11.8%, Robust: 19.5%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -860.827
## delta0: 0.3476607
## Pars: NA
## Linpars: 2.35249
#
## Fail to Reject the NULL.  No Evidence for quadratic transition
## So, if it exists, it must be linear...
##==============================================================##

#### ========  SPECIFICATION - ORDER 1 5601:6149  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 5601:6149        ##
##==============================================================##

e <- e_nab[5601:6149]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,3,0.5)
TV <- newTV(del0=1,vecpars=startpars,shape=TRshape$linear,speedopt=TRspeedopt$eta,st=ST)
# Estimate the univar series:
parScale <- c(1,1,1,1)
nDeps <- rep(1e-3,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-12)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l')

if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_1.RDS")
fName <- "TestStatDist_NAB.RDS"
# Null Hyp: Taylor Order 3
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H0,saveas=fName,numloops = 1100)
beep()
# Null Hyp: Taylor Order 1
TEST <- GenTestStatDist(e,TV,refData,testorder=TESTorder$H01,saveas=fName,numloops = 1100)
beep()

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 1 ---- obs 5601:6149         ##

## Result: Pvals=[TR2: 76.3%, Robust: 66.6%], using Taylor Order 3
## TestStat_ProbDist_NAB Logliklihood Value:  -851.28
## delta0: 0.7108208
## Pars: 
##           [,1]
## delta   2.2231890
## speed   3.2770010
## loc1    0.6011752
## loc2           NA
#
## Cannot Reject the null: H0=Model is sufficient
## So, now let's try to put them all together!
##==============================================================##


##   We have identified 4 possible transitions as follows:
# 1: e_nab[2001:3600] Lin: ## delta -1.486306    ## speed   6.999914   ## loc1    0.538707
# 2: e_nab[3601:4600] Cub: ## delta -22.8461455  ## speed   3.588742   ## loc1    0.6570495
# 3: e_nab[4601:5600] Lin: ## delta  -1.334988   ## speed   4.534029   ## loc1    0.480660
# 4: e_nab[5601:6149] Lin: ## delta   2.2231890  ## speed   3.277001   ## loc1    0.6011752

##  Approx. Locations in the full sample are:
##  0.465593, 0.6924784, 0.8264206, 0.9645544



#### ========  FINAL SPECIFICATION - ORDER 4  1:6149  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 4 ---- obs 1:6149           ##
##==============================================================##

#  NAB:
e <- e_nab
# Create Order 4 TV object with starting parameters
ST <- 1:length(e)/length(e)
startpars <- c(1,6,0.1,1,4,0.3,1,4,0.6,1,4,0.7)
Shape <- c(TRshape$linear,TRshape$cubic,TRshape$linear,TRshape$linear)
TV <- newTV(del0=1,vecpars=startpars,shape=Shape,speedopt=TRspeedopt$eta,st=ST)
# Estimate the univar series:
parScale <- c(2,2,1,2,2,1,5,2,1,5,2,1,5)
nDeps <- rep(1e-10,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-9,maxit=1000)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l')


##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 4 ---- obs 1:6149            ##
##
## Result: Pvals=[TR2: %, Robust: %], using
## TestStat_ProbDist_NAB Logliklihood Value:  -10211.48
## delta0: 4.812139
## Pars: 
##              [,1]       [,2]       [,3]        [,4]
##  delta -0.2991552 -3.7702575 20.0050056 -19.1342361
##  speed  6.9999981  4.8990590  4.6504937   4.2939839
##  loc1   0.2516205  0.3420834  0.6686268   0.7014759
##  loc2          NA         NA         NA          NA
#
## FULLY SPECIFIED!
## 
##==============================================================##


##==============================================================##
##                            THE END
##==============================================================##


