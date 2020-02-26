rm(list=ls())
gc(TRUE)

##====================== Initialisation ========================####
if (TRUE){
  #setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC  - GOOGLE DRIVE
  setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building_GW")        # Glen's home PC - GOOGLE DRIVE
  #setwd("D:/Source/Project2_Building")                                          # Anna's new Laptop
  
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
  plot(dates,e_cba,type="l",main=ptitle)
  
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
abline(h=seq(-010,010,5), v=seq(0,7000,500),col="lightgrey")
lines(e)


#### ========  SPECIFICATION - ORDER 3    ======== ####

##==============================================================##
##   TV - AUS_4Banks\CBA ---- ORDER 3  from 2018  ####
##==============================================================##

##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
##          Log the results for each test below the code.

e <- e_cba[1:5400]
e <- e_cba
# Create a TV object with delta0, pars & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,3,0.5,3,3,0.7,1,3,0.9)
startpars <- c(1,2,0.2,3,3,0.5,3,3,0.8)
TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)

# Estimate the univar series:
parScale <- c(1,1,1,0.3,3,1,0.3,3,1,0.3)
nDeps <- rep(1e-5,length(parScale))
# Estimate the univar series:
nDeps <- c(1e-5,1e-5,1e-5,1e-9,1e-5,1e-5,1e-9,1e-5,1e-5,1e-9)  #LL = -8473.838 for [1:5400]
nDeps <- c(1e-5,1e-5,1e-7,1e-9,1e-5,1e-7,1e-9,1e-5,1e-7,1e-9)  #LL = -8426.928 for [1:5400]

TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7,maxit=1000)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l',main=paste0("Loglik: ",round(TV$Estimated$value)))
## Now calculate the standard errors for the model  ##
TV <- calcStderr_TV(e,TV)
TV$Estimated$stderr

TV <- calcParamStats_TV(Tv)

TV$Estimated$parsVector
TV$Estimated$stderr
TV$Estimated$tStat
TV$Estimated$PValues

saveRDS(TV,"Estimated_TVOrd3_CBA.RDS")

TV <- readRDS("Estimated_TVOrd3_CBA.RDS")
View(TV)

if (!exists("refData")) refData <- readRDS("RefData_withGarch_CBA_2018.RDS")
fName <- "TestStatDist_CBA.RDS"
# Calculate the Reference test statistics:
refTest <- list()
refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1600)


#### ========  SPECIFICATION - ORDER 3  - Results  ======== ##


##==============================================================##
##  TV - AUS_4Banks\CBA ---- ORDER 3 ---- obs 1:7027            ##
##
## Result: Pvals=[TR2: %, Robust: %], using
## TestStat_ProbDist_CBA Logliklihood Value:  -10817.27 - Can calc StdErrs
## delta0: 1.310083
## Pars: 
##            [,1]       [,2]       [,3]
##  delta -0.6138192 8.5151909 -8.085892
##  speed  6.9999994 6.1432614  4.471730
##  loc1   0.4215485 0.5952681  0.647414
##  loc2          NA         NA         NA
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
#
#
# > TV$Estimated$parsVector
# delta0    delta1    speed1    loc1.1    delta2    speed2    loc2.1    delta3    speed3    loc3.1 
# 1.310083 -0.613819  6.999999  0.421548  8.515191  6.143261  0.595268 -8.085892  4.471730  0.647414 
# > TV$Estimated$stderr
# delta0   delta1   speed1   loc1.1   delta2   speed2   loc2.1   delta3   speed3   loc3.1 
# 0.034062 0.045020 0.237471 0.000813 0.586107 0.115458 0.000908 0.588220 0.164732 0.000215 
# > TV$Estimated$tStat
# delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
# 38.46172  -13.63436   29.47728  518.50923   14.52839   53.20775  655.58150  -13.74637   27.14549 3011.22791 
# > TV$Estimated$PValues
# delta0* delta1* speed1* loc1.1* delta2* speed2* loc2.1* delta3* speed3* loc3.1* 
#   0       0       0       0       0       0       0       0       0       0 

## I would recommend this model!
## All estimates are significant, and LL is -10817  vs. -10778 for the Ord4 Model with some non-significant params 


##==============================================================##








##==============================================================##
##  TV.ORDER 2, MTVGJR Model  ####
##==============================================================##

## e = cba[1:7027]
## delta0: 1.310083
## Pars: 
##            [,1]       [,2]       [,3]
##  delta -0.6138192 8.5151909 -8.085892
##  speed  6.9999994 6.1432614  4.471730
##  loc1   0.4215485 0.5952681  0.647414
## Logliklihood Value:  -10817.27 - Can calc StdErrs

## Start pars:
##  startpars <- c(1,3,0.5,3,3,0.7,1,3,0.9)
##  startpars <- c(1,2,0.2,3,3,0.5,3,3,0.8)

if (TRUE) {
  
  e <- e_cba
  nr.obs <- length(e)
  st <- 1:nr.obs/nr.obs
  estCtrl <- list(calcSE = TRUE,verbose = TRUE)
  
  TV <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
  
  OptCtrl <- TV$optimcontrol
  OptCtrl$reltol <- 1e-7
  #OptCtrl$ndeps <- rep(1e-5,TV@nr.pars)
  OptCtrl$ndeps <- c(1e-5, 1e-5,1e-5,1e-5, 1e-5,1e-5,1e-5, 1e-5,1e-5,1e-5)
  OptCtrl$parscale <- c(3, 1,1,0.3, 3,1,0.3, 3,1,0.3)
  TV$optimcontrol <- OptCtrl
  
  TV$delta0 <- 1
  TV$pars["deltaN",] <- c(1,3,1)
  TV$pars["speedN",] <- c(3,3,3)
  TV$pars["locN1",] <- c(0.2, 0.5, 0.8)
  
  TV <- estimateTV(e,TV,estCtrl)
  
  summary(TV)
  plot(TV)
  
}

#  Univar TV  Estimation Results:
#   
#   Delta0 = 0 se0 =  NaN 
# 
#                st1      se1 sig      st2      se2 sig.1
#   deltaN 1.2438307 0.060485 **  1.242096 0.030272   ** 
#   speedN 4.6148744 0.130642 **  3.097853 0.167118   *  
#   locN1  0.5580971 0.002868 **  0.893777 0.017572   ** 
#   locN2         NA      NaN           NA      NaN      
# 
# Log-liklihood value:  -11512.76

# Estimation Results:
#   
#   Delta0 = 1.354377 se0 =  0.03582 
# 
#               st1      se1 sig       st2      se2 sig.1        st3      se3 sig.2
#   deltaN -0.668492 0.049531 *   13.626951 3.941173   *   -13.183815 3.935306   *  
#   speedN  6.879981 0.960068 *    4.799864 0.197317   **    4.289804 0.094964   ** 
#   locN1   0.405107 0.002171 **   0.576614 0.004071   **    0.610202 0.006147   ** 
#   locN2         NA      NaN            NA      NaN               NA      NaN      
# 
# Log-liklihood value:  -11046.54

# Now Test: 
if (TRUE) {
  
  # Read in the Garch'd Reference Data
  RefData <- readRDS("Gen_Refdata/RefData_withGarch_CBA_2020.RDS")
  
  TV <- setTaylorOrder(2,TV)
  
  RefTests <- list()
  RefTests$LMTR2 <- LM.TR2(e,TV)
  RefTests$LMRobust <- LM.Robust(e,TV)
  
  simcontrol <- list()
  simcontrol$saveAs <- "CBA"
  simcontrol$numLoops <- 120
  simcontrol$numCores <- 3
  TEST <- testStatDist(TV,RefData,RefTests,simcontrol)
  
  # View Distributions
  hist(TEST$Stat_TR2,breaks = 15)
  hist(TEST$Stat_Robust,breaks = 15)  
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
  OptCtrl$reltol <- 1e-9
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
    lG1 <- MTV1$results[[1]]$garch   # For subsequent iterations - update the index!
    e <- e_orig/sqrt(lG1@h)
    
    #lT1 <- MTV1@tv                # For first iteration
    lT1 <- MTV1$results[[1]]$tv   # For subsequent iterations - update the index!
    
    summary(lT1)
    ll.target <- getTargetValue(e,lT1)
    
    lT1$delta0 <- lT1$Estimated$delta0
    lT1$pars <- lT1$Estimated$pars 
    
    lT1$pars["deltaN",] <- c(1,3,1)
    lT1$pars["speedN",] <- c(3,3,3)
    #lT1$pars["locN1",] <- c(0.1,0.4,0.5)
    lT1$pars["locN1",] <- c(0.2,0.5,0.8)
    
    OptCtrl <- TV$optimcontrol
    OptCtrl$reltol <- 1e-7
    #OptCtrl$ndeps <- rep(1e-5,TV@nr.pars)
    OptCtrl$ndeps <- c(1e-5, 1e-5,1e-5,1e-5, 1e-5,1e-5,1e-5, 1e-5,1e-5,1e-5)
    OptCtrl$parscale <- c(3, 1,1,0.3, 3,1,0.3, 3,1,0.3)
    TV$optimcontrol <- OptCtrl
    
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
    
    lG1$pars["omega",1] <- 0.02
    lG1$pars["alpha",1] <- 0.03
    lG1$pars["beta",1] <- 0.80
    lG1$pars["gamma",1] <- 0.05
    
    optCtrl <- lG1$optimcontrol
    optCtrl$reltol <- 1e-7
    optCtrl$ndeps <- c(1e-5,1e-5,1e-7,1e-5)
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
  MTV1$Estimated[[1]] <- MTV1$results[[1]]
  
  
  # Save the estimated model
  saveRDS(MTV1,"Output/CBA_mtvgjr_manual.RDS")
}

## Notes: Could not get a reliable estimate for TV in second round.

##==============================================================##
##                            THE END
##==============================================================##

MTV1 <- readRDS("Output/CBA_mtvgjr_manual.RDS")
summary(MTV1$results[[1]]$tv)
plot(MTV1$results[[1]]$tv)
summary(MTV1$results[[2]]$tv)
plot(MTV1$results[[2]]$tv)



