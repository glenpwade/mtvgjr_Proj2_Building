rm(list=ls())
gc(TRUE)

##====================== Initialisation ========================####
if (TRUE){
  #setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC & Laptop - GOOGLE DRIVE
  setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building_GW")        # Glen's home PC - GOOGLE DRIVE

  source("clsTV.r")
  library(parallel)
  library(RevoUtilsMath)
  maxThreads <- detectCores()-1
  setMKLthreads(maxThreads)
}
##====================== Initialisation ========================##


##====================== Data Load =============================####
if (TRUE){
  # Read AUS_4Banks data from saved file
  mydata <- readRDS("AUS_bankreturns_1992_2020.RDS")
  dates <- mydata$date
  
  # Create scaled-up & de-meaned individual data series vectors (Percentage Returns):
  e_wbc <- mydata$return_wbc * 100 
  e_wbc <- e_wbc - mean(e_wbc)
  e <- e_wbc
  
  # Plot the returns
  ptitle <- "WBC Returns"
  plot(dates,e,type="l",main=ptitle)
  
  # Read in the Garch'd Reference Data
  RefData <- readRDS("Gen_Refdata/RefData_withGarch_WBC_2020.RDS")
  
  # Tidy Up
  rm(mydata)
  gc(TRUE)
}
##====================== Data Load =============================##


##   
##============  General Specification Method    ===============####
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

# Plot returns using index as x-axis:
ptitle <- "WBC Returns"
plot(e,type="l",main = ptitle)
abline(h=seq(-10,10,5), v=seq(0,7000,500),col="lightgrey")
lines(e)



##==============================================================##
##  TV.ORDER 0, obs 1:7108  ####
##==============================================================##

if (TRUE) {

  e <- e_wbc
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  
  TV <- tv(st,tvshape$delta0only)
  TV <- estimateTV(e,TV)
  
  # Now Test against WBC
  TV <- setTaylorOrder(1,TV)
  
  RefTests <- list()
  RefTests$LMTR2 <- LM.TR2(e,TV)
  RefTests$LMRobust <- LM.Robust(e,TV)
  
  simcontrol <- list()
  simcontrol$saveAs <- "WBC.ORD0.1_1-7108.RDS"
  simcontrol$numLoops <- 1000
  simcontrol$numCores <- 3
  TEST <- testStatDist(TV,RefData,RefTests,simcontrol)

  # View Distributions
  hist(TEST$Stat_TR2,breaks = 15)
  hist(TEST$Stat_Robust,breaks = 15)
  
}  
##==============================================================##
##  TV - AUS_4Banks/WBC ---- ORDER 0.1 ---- obs 1:7108          ##

## Result: Pvals=[TR2: 73%,Robust= 69%] 
## Fail to Reject the null: H0=Model is sufficient
##  This is likely caused by issues with the optimiser 
##  where there are many transitions
##  So, try a smaller range... 
##==============================================================##



##==============================================================##
##  TV.ORDER 0, obs 1:4000  ####
##==============================================================##

if (TRUE) {
  
  e <- e_wbc[1:4000]
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  
  TV <- tv(st,tvshape$delta0only)
  TV <- estimateTV(e,TV)
  
  # Now Test against WBC
  TV <- setTaylorOrder(2,TV)
  
  RefTests <- list()
  RefTests$LMTR2 <- LM.TR2(e,TV)
  RefTests$LMRobust <- LM.Robust(e,TV)
  
  simcontrol <- list()
  simcontrol$saveAs <- "WBC.ORD0.2_1-4000.RDS"
  simcontrol$numLoops <- 1200
  simcontrol$numCores <- 3
  TEST <- testStatDist(TV,RefData,RefTests,simcontrol)
  
  # View Distributions
  hist(TEST$Stat_TR2,breaks = 15)
  hist(TEST$Stat_Robust,breaks = 15)
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/WBC ---- obs 1:4000          

## ORDER 0.1 
## Result: Pvals=[TR2= 3%,Robust= 4%] 
## ORDER 0.2
## Result: Pvals=[TR2= 8%,Robust= 2%] 
## ORDER 0.3
## Result: Pvals=[TR2= %,Robust= %] 

## Reject the null: H0=Model is sufficient
## Conclusion: Strong evidence of a transition here

##==============================================================##


##==============================================================##
##  TV.ORDER 1, obs 1:4000  ####
##==============================================================##

if (TRUE) {
  
  e <- e_wbc[1:4000]
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  
  # Add some estimation controls
  estCtrl <- list()
  estCtrl$calcSE <- TRUE
  estCtrl$verbose <- TRUE
  
  TV <- tv(st,tvshape$single)
  
  OptCtrl <- TV$optimcontrol
  OptCtrl$reltol <- 1e-9
  OptCtrl$ndeps <- c(1e-7,1e-7,1e-7,1e-5)
  #OptCtrl$parscale <- c(6,6,3,1)
  TV$optimcontrol <- OptCtrl
  
  TV$delta0 <- 0.5
  TV$pars[1,1] <- 0.01
  TV$pars[2,1] <- 2.05
  TV$pars[3,1] <- 0.6

  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)
  
  # Now Test against WBC
  TV <- setTaylorOrder(2,TV)
  
  RefTests <- list()
  RefTests$LMTR2 <- LM.TR2(e,TV)
  RefTests$LMRobust <- LM.Robust(e,TV)
  
  simcontrol <- list()
  simcontrol$saveAs <- "WBC.ORD1.2_1-4000.RDS"
  simcontrol$numLoops <- 1200
  simcontrol$numCores <- 3
  TEST <- testStatDist(TV,RefData,RefTests,simcontrol)
  
  # View Distributions
  hist(TEST$Stat_TR2,breaks = 15)
  hist(TEST$Stat_Robust,breaks = 15)
  
}
##==============================================================##

# Estimation Results:
#   
#   Delta0 = 0.000207 se0 =  7e-06 
# 
#             st1      se1
# deltaN -0.000113 0.000009
# speedN  2.875682 0.276979
# locN1   0.646473 0.029049
# locN2         NA      NaN
# 
# Log-liklihood value:  11811.23

##==============================================================##

##  TV - AUS_4Banks/WBC ---- ---- obs 1:5000          
## ORDER 1.1 
## Result: Pvals=[TR2= 92%,Robust= 90%] 
## ORDER 1.2
## Result: Pvals=[TR2= 40%,Robust= 43%] 
## ORDER 1.3
## Result: Pvals=[TR2= %,Robust=%] 

## Cannot Reject the null
## Conclusion: No evidence of another transition here

##==============================================================##



##==============================================================##
##  TV.ORDER 0, obs 3501:5500  ####
##==============================================================##

if (TRUE) {
  
  e <- e_wbc[3501:5500]
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  
  TV <- tv(st,tvshape$delta0only)
  
  # Add some estimation controls
  estCtrl <- list()
  estCtrl$calcSE <- TRUE
  estCtrl$verbose <- TRUE
  
  OptCtrl <- TV$optimcontrol
  OptCtrl$reltol <- 1e-11
  OptCtrl$ndeps <- OptCtrl$ndeps * 1e-2
  #OptCtrl$parscale
  TV$optimcontrol <- OptCtrl
  
  TV$delta0
  TV$pars[1,1] <- 0.05
  TV$pars[2,1] <- 2.05
  TV$pars[3,1] <- 0.65
  
  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)
  
  # Now Test against WBC
  TV <- setTaylorOrder(2,TV)
  
  RefTests <- list()
  RefTests$LMTR2 <- LM.TR2(e,TV)
  RefTests$LMRobust <- LM.Robust(e,TV)
  
  simcontrol <- list()
  simcontrol$saveAs <- "WBC.ORD0.2_3501-5500.RDS"
  simcontrol$numLoops <- 1200
  simcontrol$numCores <- 3
  TEST <- testStatDist(TV,RefData,RefTests,simcontrol)
  
  # View Distributions
  hist(TEST$Stat_TR2,breaks = 15)
  hist(TEST$Stat_Robust,breaks = 15)
  
}  # End: if (DoThis)...
##==============================================================##
##  TV - AUS_4Banks/WBC ---- ---- obs 3501:5500
## ORDER 0.1 
## Result: Pvals=[TR2= 52%,Robust= 23%] 
## ORDER 0.2
## Result: Pvals=[TR2= 3.4%,Robust= 0.4%] 
## ORDER 0.3
## Result: Pvals=[TR2= %,Robust= %] 

## Reject the null
## Conclusion: Strong evidence of a double transition here
##             Next: Lets try to estimate the 2 transitions we've found
##==============================================================##



##==============================================================##
##  TV.ORDER 2, obs 1:5500  ####
##==============================================================##

if (TRUE) {
  
  e <- e_wbc[1:5500]
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  
  # Add some estimation controls
  estCtrl <- list()
  estCtrl$calcSE <- TRUE
  estCtrl$verbose <- TRUE
  
  TV <- tv(st,c(tvshape$single,tvshape$double1loc))

  OptCtrl <- TV$optimcontrol
  OptCtrl$reltol <- 1e-9
  OptCtrl$ndeps <- c(1e-7,1e-7,1e-5,1e-5,1e-7,1e-7,1e-5)
  #OptCtrl$parscale <- c(6,6,2,1,6,2,1)
  OptCtrl$parscale <- c(3,3,1,1,3,1,1)
  TV$optimcontrol <- OptCtrl
  
  TV$delta0 <- 0.1
  TV$pars["deltaN",1] <- 0.05
  TV$pars["speedN",1] <- 3.05
  TV$pars["locN1",1] <- 0.3
  TV$pars["deltaN",2] <- 0.05
  TV$pars["speedN",2] <- 3.05
  TV$pars["locN1",2] <- 0.6
  
  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)
  
  # Now Test against WBC
  TV <- setTaylorOrder(2,TV)
  
  RefTests <- list()
  RefTests$LMTR2 <- LM.TR2(e,TV)
  RefTests$LMRobust <- LM.Robust(e,TV)
  
  simcontrol <- list()
  simcontrol$saveAs <- "WBC.ORD2sd.2_1-5500.RDS"
  simcontrol$numLoops <- 1200
  simcontrol$numCores <- 3
  TEST <- testStatDist(TV,RefData,RefTests,simcontrol)
  
  # View Distributions
  hist(TEST$Stat_TR2,breaks = 15)
  hist(TEST$Stat_Robust,breaks = 15)
  
}

##==============================================================##

# Estimation Results:
#   
#   Delta0 = 0.001825 se0 =  0.000163 
# 
#             st1       st2      se1      se2
# deltaN -0.000072 -0.001614 0.000008 0.000163
# speedN  3.461076  6.999988 0.418388 0.175481
# locN1   0.397173  0.774516 0.022339 0.003241
# locN2         NA        NA      NaN      NaN
# 
# Log-liklihood value:  15732.82 

##==============================================================##

##==============================================================##
##  TV - AUS_4Banks/WBC ---- ---- obs 1:5500          
## ORDER 2.1 
## Result: Pvals=[TR2= 32%,Robust= 16%] 
## ORDER 2.2
## Result: Pvals=[TR2= 39%,Robust= 33%] 

## Fail to Reject the null
## Conclusion: No Evidence of another transition here

##==============================================================##



##==============================================================##
##  TV.ORDER 2, obs 1:7108  ####
##==============================================================##

if (TRUE) {
  
  e <- e_wbc[1:7108]
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  
  # Add some estimation controls
  estCtrl <- list()
  estCtrl$calcSE <- TRUE
  estCtrl$verbose <- TRUE
  
  TV <- tv(st,c(tvshape$single,tvshape$double1loc))
  
  OptCtrl <- TV$optimcontrol
  OptCtrl$reltol <- 1e-9
  OptCtrl$ndeps <- c(1e-7,1e-7,1e-5,1e-5,1e-7,1e-7,1e-5)
  OptCtrl$parscale <- c(6,6,1,1,6,1,1)
  #OptCtrl$parscale <- c(3,3,1,1,3,1,1)
  TV$optimcontrol <- OptCtrl
  
  TV$delta0 <- 0.1
  TV$pars["deltaN",1] <- 0.05
  TV$pars["speedN",1] <- 3.05
  TV$pars["locN1",1] <- 0.3
  TV$pars["deltaN",2] <- 0.05
  TV$pars["speedN",2] <- 3.05
  TV$pars["locN1",2] <- 0.6
  
  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)
  
  # Now Test against WBC
  TV <- setTaylorOrder(2,TV)
  
  RefTests <- list()
  RefTests$LMTR2 <- LM.TR2(e,TV)
  RefTests$LMRobust <- LM.Robust(e,TV)
  
  simcontrol <- list()
  simcontrol$saveAs <- "WBC.ORD2sd.2_1-7108.RDS"
  simcontrol$numLoops <- 1200
  simcontrol$numCores <- 4
  TEST <- testStatDist(TV,RefData,RefTests,simcontrol)
  
  # View Distributions
  hist(TEST$Stat_TR2,breaks = 15)
  hist(TEST$Stat_Robust,breaks = 15)
  
}

##==============================================================##

# Estimation Results:
#   
#   Delta0 = 0.001618 se0 =  0.000105 
# 
#             st1       st2      se1      se2
# deltaN -0.000079 -0.001407 0.000008 0.000104
# speedN  2.963292  7.000000 0.252549 0.000000
# locN1   0.319891  0.606221 0.019991 0.001652
# locN2         NA        NA      NaN      NaN
# 
# Log-liklihood value:  20590.96

##==============================================================##

##==============================================================##
##  TV - AUS_4Banks/WBC ---- ---- obs 1:7108          
## ORDER 2.1 
## Result: Pvals=[TR2= 47%,Robust= 93%] 
## ORDER 2.2
## Result: Pvals=[TR2= NaN%,Robust= 99.8%] 

## Fail to Reject the null
## Conclusion: No Evidence of another transition here
##             Model is fully specified
##==============================================================##


##==============================================================##
##  TV.ORDER 2, MTVGJR Model  ####
##==============================================================##

##====================== Data Load =============================####
if (TRUE){
  # Read AUS_4Banks data from saved file
  mydata <- readRDS("AUS_bankreturns_1992_2020.RDS")
  dates <- mydata$date
  
  # Create scaled-up & de-meaned individual data series vectors (Percentage Returns):
  e_wbc <- mydata$return_wbc * 100 
  e_wbc <- e_wbc - mean(e_wbc)
  e <- e_wbc
  
  # Plot the returns
  ptitle <- "WBC Returns"
  plot(dates,e,type="l",main=ptitle)
  
  # Read in the Garch'd Reference Data
  RefData <- readRDS("Gen_Refdata/RefData_withGarch_WBC_2020.RDS")
  
  # Tidy Up
  rm(mydata)
  gc(TRUE)
}

if (TRUE) {
  
  source("clsMTVGARCH.r")
  
  e <- e_wbc
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  
  TV <- tv(st,c(tvshape$single,tvshape$double1loc))
  
  OptCtrl <- TV$optimcontrol
  OptCtrl$reltol <- 1e-9
  OptCtrl$ndeps <- c(1e-7,1e-7,1e-5,1e-5,1e-7,1e-7,1e-5)
  OptCtrl$parscale <- c(6,6,1,1,6,1,1)
  TV$optimcontrol <- OptCtrl
  
  TV$delta0 <- 1
  TV$pars["deltaN",1] <- 1
  TV$pars["speedN",1] <- 3.05
  TV$pars["locN1",1] <- 0.3
  TV$pars["deltaN",2] <- 1
  TV$pars["speedN",2] <- 3.05
  TV$pars["locN1",2] <- 0.6
  
  estCtrl <- list(calcSE = TRUE,verbose = TRUE)
  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)

  # Create the MTV-GJR object
  MTV_wbc <- mtvgarch(TV,garchType = garchtype$gjr)
  
  # Estimate the MTV-GJR object
  MTV_wbc <- estimateMTVGARCH(e,MTV_wbc)  
  
  # Look at the univariate MTV-GJR object:
  # Plot univar TV 'g'
  plot(MTV_wbc$tv)
  # Overlay the final estimated TV 'g'
  lines(MTV_wbc$Estimated$tv@g,type="l",col="red")
  # If the scale is off, plot the Estimated TV 'g'
  plot(MTV_wbc$Estimated$tv@g,type="l",col="red")
  
  # Plot univar Garch 'h'
  plot(MTV_wbc$garch@h,type="l")
  # Overlay the final estimated Garch 'h'
  lines(MTV_wbc$Estimated$garch@h,type="l",col="red")
  # If the scale is off, plot the Estimated Garch 'h'
  plot(MTV_wbc$Estimated$garch@h,type="l",col="red")
  
  # Look at the final parameters
  summary(MTV_wbc$Estimated$tv)
  summary(MTV_wbc$Estimated$garch)
  
  # Save the estimated model
  saveRDS(MTV_wbc,"Output/WBC_mtvgjr.RDS")
  

}