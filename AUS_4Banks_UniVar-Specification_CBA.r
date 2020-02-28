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
  e_cba <- mydata$return_cba * 100 
  e_cba <- e_cba - mean(e_cba)
  e <- e_cba
  
  # Plot the returns
  ptitle <- "CBA Returns"
  plot(dates,e,type="l",main=ptitle)
  
  # Read in the Garch'd Reference Data
  RefData <- readRDS("Gen_Refdata/RefData_withGarch_CBA_2020.RDS")
  
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
ptitle <- "CBA Returns"
plot(e,type="l",main = ptitle)
abline(h=seq(-10,10,5), v=seq(0,7000,500),col="lightgrey")
lines(e)




##==============================================================##
##  TV.ORDER 2, obs 1:7108  ####
##==============================================================##

if (TRUE) {
  
  e <- e_cba[1:7108]
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  
  # Add some estimation controls
  estCtrl <- list()
  estCtrl$calcSE <- TRUE
  estCtrl$verbose <- TRUE
  
  TV <- tv(st,c(tvshape$single,tvshape$double1loc))
  
  OptCtrl <- TV$optimcontrol
  OptCtrl$reltol <- 1e-9
  OptCtrl$ndeps <- c(1e-7,1e-7,1e-7,1e-7,1e-7,1e-7,1e-7)
  #OptCtrl$parscale <- c(6,3,1,1,3,1,1)
  OptCtrl$parscale <- c(3,3,1,1,3,1,1)
  TV$optimcontrol <- OptCtrl
  
  TV$delta0 <- 1
  TV$pars["deltaN",] <- c(1,1)
  TV$pars["speedN",1] <- 3.05
  TV$pars["locN1",1] <- 0.2
  
  TV$pars["speedN",2] <- 3.05
  TV$pars["locN1",2] <- 0.6

  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)
  
  # Fix the std errors
  h0 <- TV$Estimated$hessian
  h1 <- h0[-6,-6]
  se1 <- try(stdErrors <- sqrt(-diag(qr.solve(h1))))
  h2 <- h1[-3,-3]
  se2 <- try(stdErrors <- sqrt(-diag(qr.solve(h2))))
  
  
  # Now Test against CBA
  TV <- setTaylorOrder(1,TV)
  
  RefTests <- list()
  RefTests$LMTR2 <- LM.TR2(e,TV)
  RefTests$LMRobust <- LM.Robust(e,TV)
  
  simcontrol <- list()
  simcontrol$saveAs <- "CBA.ORD2.2_1-7108.RDS"
  simcontrol$numLoops <- 120
  simcontrol$numCores <- 3
  TEST <- testStatDist(TV,RefData,RefTests,simcontrol)
  
  # View Distributions
  hist(TEST$Stat_TR2,breaks = 15)
  hist(TEST$Stat_Robust,breaks = 15)
  
}

##==============================================================##

# Estimation Results:
#   
#   Delta0 = 15.2539 se0 =  0.061529 
# 
#               st1      se1 sig        st2      se2   sig.1
#   deltaN 0.433884 0.051305 *   -14.189119 0.067410   ***
#   speedN 6.999964 0.182940 **    7.000000 0.000000   ***
#   locN1  0.201900 0.001770 **    0.601815 0.001578   ***
#   locN2        NA      NaN             NA      NaN      
# 
# Log-liklihood value:  -11794.82
##==============================================================##


##==============================================================##
##  TV - AUS_4Banks/CBA ---- ---- obs 1:7108          
## ORDER 2.1 
## Result: Pvals=[TR2= 55%,Robust= 26%] 
## ORDER 2.2
## Result: Pvals=[TR2= 62%,Robust= 51%] 
## ORDER 2.3
## Result: Pvals=[TR2= %,Robust= %] 

## Reject the null
## Conclusion: No Evidence of another transition here
##             Model is fully specified
##==============================================================##


## From 2018:
##==============================================================##
##  TV - AUS_4Banks\CBA ---- ORDER 3 ---- obs 1:7027            ##
##
## Result: Pvals=[TR2: 78.5%, Robust: 54.9%], using
## TestStat_ProbDist_CBA Logliklihood Value:  -11438.12
## delta0: 2.263294
## Pars: 
##            [,1]       [,2]       [,3]
##  delta -1.3629961 23.1683258 -22.6972043
##  speed  6.0774728  4.7347456   4.2034498
##  loc1   0.1957340  0.6023369   0.6273989
##  loc2   0.4019763         NA          NA
##
##==============================================================##

if(TRUE){
  e <- e_cba
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  estCtrl <- list(calcSE = TRUE,verbose = TRUE)
  
  TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
  
  OptCtrl$fnscale <- -1
  OptCtrl$reltol <- 1e-9
  OptCtrl$parscale <- c(3, 3,1,1, 3,3,1, 3,3,1)
  OptCtrl$ndeps <- rep(1e-5, TV@nr.pars)
  #OptCtrl$ndeps <- rep(1e-5,TV@nr.pars)
  #OptCtrl$ndeps <- c(1e-5, 1e-5,1e-7,1e-5, 1e-7,1e-5,1e-5, 1e-7,1e-5,1e-5)
  TV$optimcontrol <- OptCtrl
  
  TV$delta0 <- 1.0
  TV$pars["deltaN",] <- c(1, 3, -3)
  TV$pars["speedN",] <- c(1, 3, 3)
  TV$pars["locN1",] <- c(0.3, 0.6, 0.7)
  #TV$pars["locN2",] <- c(0.3, NA, NA)
  
  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)

}

# Estimation Results:
#   
#   Delta0 = 1.349635 se0 =  0.035672 
# 
# st1      se1 sig       st2      se2 sig.1        st3      se3 sig.2
# deltaN -0.667736 0.048937 *   13.362309 3.673008   *   -12.917186 3.667844   *  
#   speedN  6.999993 0.000000 ***  4.823184 0.194542   **    4.290551 0.095564   ** 
#   locN1   0.404800 0.001835 ***  0.576191 0.003887   **    0.610496 0.006018   ** 
#   locN2         NA      NaN            NA      NaN               NA      NaN      
# 
# Log-liklihood value:  -11035.09




##==============================================================##
##  TV.ORDER 3, MTVGJR Model - Manual Estimation   ####
##==============================================================##

if (TRUE) {
  
  # Create e, TV1 & G1 
  # make an original copy of the data:
  e_orig <- e_cba
  
  estCtrl <- list(calcSE = TRUE, verbose = TRUE)
  
  TV1 <- TV  # Estimated above
  e <- e_orig/sqrt(TV1@g)
  
  G1 <- garch(garchtype$gjr)
  
  OptCtrl <- G1$optimcontrol
  OptCtrl$reltol <- 1e-7
  OptCtrl$ndeps <- c(1e-5,1e-5,1e-5,1e-5)
  OptCtrl$parscale <- c(1,1,9,1)
  G1$optimcontrol <- OptCtrl
  
  G1$pars["omega",1] <- 0.02
  G1$pars["alpha",1] <- 0.04
  G1$pars["beta",1] <- 0.8
  G1$pars["gamma",1] <- 0.03
  
  G1 <- estimateGARCH(e,G1,estCtrl)
  summary(G1)
  plot(G1)
  
  MTV1 <- mtvgarch(e,TV1,G1)
  
  ###  START THE PAIR_WISE ESTIMATION HERE  ###
  
  # Filter out Garch & estimate TV component
  if(TRUE) {
    
    lG1 <- MTV1@garch                # For first iteration
    #lastResult <- length(MTV1$results)
    #lG1 <- MTV1$results[[lastResult]]$garch   # For subsequent iterations 
    e <- e_orig/sqrt(lG1@h)
    
    lT1 <- MTV1@tv                # For first iteration
    #lT1 <- MTV1$results[[lastResult]]$tv   # For subsequent iterations 
    
    summary(lT1)
    plot(lT1)
    ll.target <- getTargetValue(e,lT1)
    
    lT1$delta0 <- lT1$Estimated$delta0
    lT1$pars <- lT1$Estimated$pars 
    
    lT1$pars["deltaN",] <- c(-1, 12, -11)
    lT1$pars["speedN",] <- c(5.5, 3, 3)
    lT1$pars["locN1",] <- c(0.3, 0.4, 0.5)
    
    OptCtrl <- lT1$optimcontrol
    OptCtrl$reltol <- 1e-7
    #OptCtrl$parscale <- c(3,3,0.05, 5,3,0.05, 5,3,0.05)
    OptCtrl$parscale <- c(1,3,1, 9,3,1, 9,3,1)
    
    #OptCtrl$ndeps <- rep(1e-5,lT1@nr.pars)
    OptCtrl$ndeps <- c(1e-7,1e-7,1e-7, 1e-5,1e-7,1e-7, 1e-5,1e-7,1e-7)
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
    lG1$pars <- lG1$Estimated$pars
    
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
  MTV1$Estimated[[1]] <- MTV1$results[[2]]
  
  # Save the estimated model
  saveRDS(MTV1,"Output/CBA_mtvgjr_manual.RDS")
}







##==============================================================##
##                            THE END
##==============================================================##



# lT1 <- MTV1@tv                # For first iteration
# #lT1 <- MTV1$results[[lastResult]]$tv   # For subsequent iterations 
# 
# summary(lT1)
# plot(lT1)
# ll.target <- getTargetValue(e,lT1)
# 
# lT1$delta0 <- lT1$Estimated$delta0
# lT1$pars <- lT1$Estimated$pars 
# 
# lT1$pars["deltaN",] <- c(-1, 12, -11)
# lT1$pars["speedN",] <- c(5.5, 3, 3)
# lT1$pars["locN1",] <- c(0.3, 0.4, 0.5)
# 
# OptCtrl <- lT1$optimcontrol
# OptCtrl$reltol <- 1e-7
# #OptCtrl$parscale <- c(3,3,0.05, 5,3,0.05, 5,3,0.05)
# OptCtrl$parscale <- c(1,3,1, 9,3,1, 9,3,1)
# 
# #OptCtrl$ndeps <- rep(1e-5,lT1@nr.pars)
# OptCtrl$ndeps <- c(1e-7,1e-7,1e-7, 1e-5,1e-7,1e-7, 1e-5,1e-7,1e-7)
# lT1$optimcontrol <- OptCtrl

# Estimation Results:
#   
#   Delta0 = 1.349635 se0 =  NaN 
# 
# st1      se1 sig       st2      se2 sig.1        st3      se3 sig.2
# deltaN -0.636138 0.034757 *   22.176415 4.649566   *   -21.793676 4.655014   *  
#   speedN  6.774782 0.195077 **   4.578019 0.029249   **    4.309796 0.091418   ** 
#   locN1   0.405353 0.002594 **   0.584026 0.002503   ***   0.604268 0.003920   ** 
#   locN2         NA      NaN            NA      NaN               NA      NaN      
# 
# Log-liklihood value:  -11033.09
