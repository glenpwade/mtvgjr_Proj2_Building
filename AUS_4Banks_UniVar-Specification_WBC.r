rm(list=ls())
gc(TRUE)

##====================== Initialisation ========================####
if (TRUE){
  #setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC & Laptop - GOOGLE DRIVE
  #setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building_GW")     # Glen's home PC - GOOGLE DRIVE
  #setwd("D:/Source/Project2_Building")                                 # Anna's new Laptop
  setwd("C:/Source/MTVGJR_MGARCH/Project2_Building")
  
  source("clsTV.r")
  source("clsGARCH.r")
  source("clsMTVGARCH.r")
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



TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-6)

#### ========  SPECIFICATION - ORDER 3  - Results  ======== ####

##======================= 2018 =================================##
##  TV - AUS_4Banks\WBC ---- ORDER 3 ---- obs 1:7027            ##
##
## Result: Pvals=[TR2: %, Robust: %], using
## TestStat_ProbDist_WBC Logliklihood Value:  -11831.69
## delta0: 2.315575
## Pars: 
##            [,1]       [,2]       [,3]
##  delta -3.4379046 16.6098169 -14.4907300
##  speed  1.7795618  4.6358355   5.4461746
##  loc1   0.5571464  0.6074963   0.6336398
##  loc2          NA         NA          NA
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
#
#
# > TV$Estimated$parsVector
# delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
# 2.315575  -3.437905   1.779562   0.557146  16.609817   4.635835   0.607496 -14.490730   5.446175   0.633640 
# > TV$Estimated$stderr
# delta0   delta1   speed1   loc1.1   delta2   speed2   loc2.1   delta3   speed3   loc3.1 
# 0.171320 0.344231 0.175136 0.028286 3.901832 0.115704 0.004249 3.912726 0.399028 0.002118 
# > TV$Estimated$tStat
# delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
# 13.516080  -9.987203  10.161029  19.696882   4.256928  40.066333 142.973876  -3.703487  13.648604 299.169027 
# > TV$Estimated$PValues
# delta0*  delta1*  speed1*  loc1.1*  delta2*  speed2*  loc2.1*  delta3*  speed3*  loc3.1* 
#   0.000000 0.000000 0.000000 0.000000 0.000020 0.000000 0.000000 0.000204 0.000000 0.000000 
# 
##  I would recommend using this model!
#
##==============================================================##


if (TRUE) {
  
  e <- e_wbc
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  estCtrl <- list(calcSE = TRUE,verbose = TRUE)
  
  TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
  
  OptCtrl <- TV$optimcontrol
  OptCtrl$reltol <- 1e-7
  OptCtrl$parscale <- c(1, 1,1,0.2, 3,1,0.2, 3,1,0.2)
  OptCtrl$ndeps <- rep(1e-5,TV@nr.pars)
  #OptCtrl$ndeps <- c(1e-5, 1e-5,1e-7,1e-5, 1e-7,1e-5,1e-5, 1e-7,1e-5,1e-5)
  TV$optimcontrol <- OptCtrl
  
  TV$delta0 <- 1.0
  TV$pars["deltaN",] <- c(1, 3, -3)
  TV$pars["speedN",] <- c(3, 3, 3)
  TV$pars["locN1",] <- c(0.3, 0.4, 0.8)

  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)

  # Save the estimated model
  saveRDS(MTV_wbc,"Output/WBC_mtvgjr.RDS")
  

}


##==============================================================##
##  TV.ORDER 3, MTVGJR Model - Manual Estimation   ####
##==============================================================##

if (TRUE) {
  
  # Create e, TV1 & G1 
  # make an original copy of the data:
  e_orig <- e
  estCtrl <- list(calcSE = TRUE, verbose = TRUE)
  
  TV1 <- TV  # Estimated above
  
  G1 <- garch(garchtype$gjr)
  
  OptCtrl <- G1$optimcontrol
  OptCtrl$reltol <- 1e-7
  OptCtrl$ndeps <- c(1e-5,1e-5,1e-5,1e-5)
  OptCtrl$parscale <- c(1,1,9,1)
  G1$optimcontrol <- OptCtrl
  
  G1$pars["omega",1] <- 0.01
  G1$pars["alpha",1] <- 0.02
  G1$pars["beta",1] <- 0.8
  G1$pars["gamma",1] <- 0.03
  
  G1 <- estimateGARCH(e,G1,estCtrl)
  summary(G1)
  plot(G1)
  
  MTV1 <- mtvgarch(e,TV1,G1)
  
  ###  START THE PAIR_WISE ESTIMATION HERE  ###
  
  # Filter out Garch & estimate TV component
  if(TRUE) {
    
    #lG1 <- MTV1@garch                # For first iteration

    lastResult <- length(MTV1$results)
    lG1 <- MTV1$results[[lastResult]]$garch   # For subsequent iterations 

    e <- e_orig/sqrt(lG1@h)
    
    #lT1 <- MTV1@tv                # For first iteration
    lT1 <- MTV1$results[[lastResult]]$tv   # For subsequent iterations 
    
    summary(lT1)
    plot(lT1)
    ll.target <- getTargetValue(e,lT1)
    
    lT1$delta0 <- lT1$Estimated$delta0
    lT1$pars <- lT1$Estimated$pars * 0.9943
    
    # lT1$pars["deltaN",] <- c(1, 1, 1)
    lT1$pars["speedN",] <- c(0.8, 3, 3)
    # lT1$pars["locN1",] <- c(0.2, 0.5, 0.5)
    
    OptCtrl <- lT1$optimcontrol
    OptCtrl$reltol <- 1e-9
    OptCtrl$parscale <- c(1,1,0.2, 3,2,0.2, 3,2,0.2)
    #OptCtrl$ndeps <- rep(1e-5,lT1@nr.pars)
    OptCtrl$parscale <- c(1,2,1, 3,3,1, 3,3,1)
    OptCtrl$ndeps <- c(1e-7,1e-7,1e-7, 1e-7,1e-7,1e-7, 1e-7,1e-7,1e-7)
    lT1$optimcontrol <- OptCtrl
    
    lT1 <- estimateTV(e,lT1,estCtrl)
    summary(lT1)  
    plot(lT1)
    
    # Check that the LL is higher than the target:
    lT1$Estimated$value > ll.target 
    
  }
  
  # Filter out TV & estimate GARCH component
  if(TRUE) {
    
    e <- e_orig/sqrt(lT1@g)
    
    summary(lG1)
    ll.target <- getTargetValue(e,garchObj = lG1)
    lG1$pars <- lG1$Estimated$pars * 0.987
    
    # lG1$pars["omega",1] <- 0.05
    # lG1$pars["alpha",1] <- 0.05
    # lG1$pars["beta",1] <- 0.80
    # lG1$pars["gamma",1] <- 0.05
    
    optCtrl <- lG1$optimcontrol
    optCtrl$reltol <- 1e-9
    optCtrl$ndeps <- c(1e-7,1e-7,1e-7,1e-7)
    optCtrl$parscale <- c(1,1,9,1)
    lG1$optimcontrol <- optCtrl
    
    lG1 <- estimateGARCH(e,lG1,estCtrl)
    summary(lG1)  
    plot(lG1)
    
    # Check that the LL is higher than the target:
    lG1$Estimated$value > ll.target 
    
  }
  
  # Now call the calculate liklihood method for mtvgarch
  
  MTV1 <- calcLoglik(MTV1,lT1,lG1)
  
  # Pick the best specification & write it to the Estimated List:
  MTV1$Estimated[[1]] <- MTV1$results[[7]]
  
  # Save the estimated model
  saveRDS(MTV1,"Output/WBC_mtvgjr_manual.RDS")
}

##==============================================================##
##                            THE END
##==============================================================##








