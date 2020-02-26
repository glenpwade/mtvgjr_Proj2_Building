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
  e_anz <- mydata$return_anz * 100 
  e_anz <- e_anz - mean(e_anz)
  e <- e_anz
  
  # Plot the returns
  ptitle <- "ANZ Returns"
  plot(dates,e_anz,type="l",main=ptitle)
  
  # Read in the Garch'd Reference Data
  RefData <- readRDS("Gen_Refdata/RefData_withGarch_ANZ_2020.RDS")
  
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
ptitle <- "ANZ Returns"
plot(e_anz,type="l",main = ptitle)
abline(h=seq(-10,10,5), v=seq(0,7000,500),col="lightgrey")
lines(e_anz)



#### ========  SPECIFICATION - ORDER 3    ======== ####

##==============================================================##
##   TV - AUS_4Banks\ANZ ---- ORDER 3 
##==============================================================##

##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
##          Log the results for each test below the code.

e <- e_anz
# Create a TV object with just delta0, pars & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,3,0.4,5,3,0.5,-5,3,0.7)
TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)

# Estimate the univar series:
startpars <- c(1,3,0.4,5,3,0.5,-5,3,0.7)
parScale <- c(1,1,1,0.3,3,1,0.3,3,1,0.3)
nDeps <- rep(1e-5,length(parScale))
#nDeps <- c(1e-5,1e-5,1e-7,1e-9,1e-5,1e-7,1e-9,1e-5,1e-7,1e-9)  

TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l',main=paste0("Loglik: ",round(TV$Estimated$value)))

## Now calculate the standard errors for the model  ##
TV <- calcStderr_TV(e,TV)
TV$Estimated$stderr

TV <- calcParamStats_TV(TV)

TV$Estimated$parsVector
TV$Estimated$stderr
TV$Estimated$tStat
TV$Estimated$PValues

saveRDS(TV,"Estimated_TVOrd3_ANZ.RDS")

TV <- readRDS("Estimated_TVOrd3_ANZ.RDS")
View(TV)


if (!exists("refData")) refData <- readRDS("RefData_withGarch_ANZ_2018.RDS")
fName <- "TestStatDist_ANZ.RDS"
# Calculate the Reference test statistics:
refTest <- list()
refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1600)

# Confirm we managed to execute more than 1000 samples for our distribution:
TEST$FailCount_TR2
TEST$FailCount_Robust



#### ========  SPECIFICATION - ORDER 3  - Results  ======== ####

##==============================================================##
##  TV - AUS_4Banks\ANZ ---- ORDER 3 ---- obs 1:7027            ##
##
## Result: Pvals=[TR2: %, Robust: %], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -11972.38
## delta0: 2.203641
## Pars: 
##            [,1]       [,2]       [,3]
##  delta -1.4643366 28.9699230 -28.291274
##  speed  6.9999574  4.7740895   4.393625
##  loc1   0.4191651  0.6068634   0.626848
##  loc2          NA         NA         NA
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
##
##
# TV$Estimated$parsVector
# delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
# 2.207522  -1.464337   6.999957   0.419165  28.969923   4.774089   0.606863 -28.291274   4.393625   0.626848 
#
# TV$Estimated$stderr
# delta0    delta1    speed1    loc1.1    delta2    speed2    loc2.1    delta3    speed3    loc3.1 
# 0.057713  0.067530  0.674091  0.001124 29.240851  0.324743  0.008400 29.233087  0.104863  0.012158 
#
# TV$Estimated$tStat
# delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
# 38.249994 -21.684244  10.384291 372.922598   0.990735  14.701130  72.245595  -0.967783  41.898715  51.558480 
#
# TV$Estimated$PValues
# delta0*  delta1*  speed1*  loc1.1*   delta2  speed2*  loc2.1*   delta3  speed3*  loc3.1* 
#   0.000000 0.000000 0.000000 0.000000 0.305757 0.000000 0.000000 0.316527 0.000000 0.000000 
##==============================================================##





##  Just for fun we'll look at a slightly simpler model, because the 2nd & 3rd transition are so similar in size & speed

##==============================================================##
##   TV - AUS_4Banks\ANZ ---- ORDER 2 ---- Single + double 
##==============================================================##

e <- e_anz
# Create a TV object with just delta0, pars & st:
ST <- 1:length(e)/length(e)
startpars <- c(1,3,0.3,1,3,0.5,0.8)
TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$double),speedopt=TRspeedopt$eta)
# Estimate the univar series:
parScale <- c(2,2,2,0.1,2,2,0.1,0.1)
nDeps <- rep(5e-5,length(parScale))
TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-5)
TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
# Estimation Outputs:
TV$Estimated$value
TV$Estimated$delta0
TV$Estimated$pars
plot(TV$Estimated$condvars,type = 'l',main=paste0("Loglik: ",round(TV$Estimated$value)))


## Now calculate the standard errors for the model  ##
TV <- calcStderr_TV(e,TV)
TV$Estimated$stderr

TV <- calcParamStats_TV(TV)

TV$Estimated$parsVector
TV$Estimated$stderr
TV$Estimated$tStat
TV$Estimated$PValues

saveRDS(TV,"Estimated_TVOrd2_ANZ.RDS")

TV <- readRDS("Estimated_TVOrd2_ANZ.RDS")
View(TV)


if (!exists("refData")) refData <- readRDS("RefData_withGarch_ANZ_2018.RDS")
fName <- "TestStatDist_ANZ.RDS"
# Calculate the Reference test statistics:
refTest <- list()
refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1600)


#### ========  SPECIFICATION - ORDER 2 (single+double)  - Results  ======== ####

##==============================================================##
##  TV - AUS_4Banks\ANZ ---- ORDER 2 (single,double) ---- obs 1:7027            ##
##
## Result: Pvals=[TR2: %, Robust: %], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -12074.9  -  Can't calculate StdErrs
## delta0: 19.29552
## Pars: 
##            [,1]       [,2]       
##  delta -1.1376929 -16.9684635
##  speed  3.6569726   6.9999143
##  loc1   0.3579471   0.6245663
##  loc2          NA   0.6351790
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
##==============================================================##

##==============================================================##
##  TV - AUS_4Banks\ANZ ---- ORDER 2 (single,double) ---- obs 1:7027            ##
##
## Result: Pvals=[TR2: %, Robust: %], using
## TestStat_ProbDist_ANZ Logliklihood Value:  -12089.84  -  Can calculate StdErrs!
## delta0: 14.89402
## Pars: 
##            [,1]       [,2]       
##  delta -0.9781297 -12.6585744
##  speed  3.3113084   6.9985303
##  loc1   0.3248167   0.6120492
##  loc2          NA   0.6551738
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
##
#
# TV$Estimated$parsVector
# delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     loc2.2
# 14.894015  -0.978130   3.311308   0.324817 -12.658574   6.998530   0.612049   0.655174
#
# > TV$Estimated$stderr
# delta0   delta1   speed1   loc1.1   delta2   speed2   loc2.1   loc2.2
# 1.846493 0.075077 0.157442 0.013135 1.845271 0.090220 0.005091 0.006079
#
# > TV$Estimated$tStat
# delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     loc2.2
# 8.066110 -13.028358  21.031923  24.729121  -6.860008  77.571824 120.221764 107.776608
#
# > TV$Estimated$PValues
# delta0* delta1* speed1* loc1.1* delta2* speed2* loc2.1* loc2.2*
#   0       0       0       0       0       0       0       0


#### ======================  Conclusion  ======================= ####

## Recommend the Order 3 model



##==============================================================##
##  MTVGJR Model - Manual Estimation   ####
##==============================================================##

if (TRUE) {
  
  ##==============================================================##
  ##  TV - AUS_4Banks\ANZ ---- ORDER 3 ---- obs 1:7027            ##
  ##
  ## Result: Pvals=[TR2: %, Robust: %], using
  ## TestStat_ProbDist_ANZ Logliklihood Value:  -11972.38
  ## delta0: 2.203641
  ## Pars: 
  ##            [,1]       [,2]       [,3]
  ##  delta -1.4643366 28.9699230 -28.291274
  ##  speed  6.9999574  4.7740895   4.393625
  ##  loc1   0.4191651  0.6068634   0.626848
  ##  loc2          NA         NA         NA
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  
  # Create e, TV1 & G1 
  # make an original copy of the data:
  e_orig <- e
  estCtrl <- list(calcSE = TRUE, verbose = TRUE)
  st <- 1:length(e)/length(e)
  
  TV1 <- tv(st,c(tvshape$single,tvshape$single,tvshape$single))
  #startpars <- c(1,3,0.4,5,3,0.5,-5,3,0.7)
  #parScale <- c(1,1,1,0.3,3,1,0.3,3,1,0.3)
  
  TV1$delta0 <- 1
  TV1$pars["deltaN",] <- c(1.0,5.0,-5.0)
  TV1$pars["speedN",] <- c(3.0,3.0,3.0)
  TV1$pars["locN1",] <- c(0.4,0.5,0.7)
  
  optCtrl <- TV1$optimcontrol
  optCtrl$reltol <- 1e-9
  optCtrl$ndeps <- c(1e-5, 1e-7,1e-7,1e-5, 1e-5,1e-5,1e-5, 1e-5,1e-5,1e-5)
  optCtrl$parscale <- c(1, 2,1,1, 3,2,1, 3,1,1)
  TV1$optimcontrol <- optCtrl
  
  TV1 <- estimateTV(e,TV1,estCtrl)
  summary(TV1)  
  plot(TV1)
  
  
  # TV OBJECT
  # 
  # Estimation Results:
  #   
  # Delta0 = 1.437446 se0 =  0.182108 
  # 
  #               st1      se1 sig       st2      se2   sig.1      st3      se3   sig.2
  #   deltaN 0.449112 0.187388 *   12.204470 1.240690   *   -12.679960 1.242423   *  
  #   speedN 6.929265 0.148397 **   6.999988 0.000000   ***   4.189196 0.105238   ** 
  #   locN1  0.018103 0.002142 *    0.572519 0.000452   ***   0.612876 0.002874   ***
  #   locN2        NA      NaN            NA      NaN               NA      NaN      
  # 
  # Log-liklihood value:  -12361.55
  
  
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
    
    lG1 <- MTV1@garch                # For first iteration
    lG1 <- MTV1$results[[1]]$garch   # For subsequent iterations - update the index!
    e <- e_orig/sqrt(lG1@h)
    
    #lT1 <- MTV1@tv                # For first iteration
    lT1 <- MTV1$results[[1]]$tv   # For subsequent iterations - update the index!
    
    summary(lT1)
    ll.target <- getTargetValue(e,lT1)
    
    lT1$delta0 <- lT1$Estimated$delta0
    lT1$pars <- lT1$Estimated$pars 
    
    lT1$pars["deltaN",] <- c(1.0,1.0)
    lT1$pars["speedN",] <- c(3.0,3.0)
    lT1$pars["locN1",] <- c(0.2,0.4)
    #lT1$pars["locN2",] <- 0
    
    optCtrl <- lT1$optimcontrol
    optCtrl$reltol <- 1e-9
    optCtrl$ndeps <- c(1e-7,1e-5,1e-7,1e-7,1e-5,1e-7)
    optCtrl$parscale <- c(3,1,1,3,1,1)
    lT1$optimcontrol <- optCtrl
    
    lT1 <- estimateTV(e,lT1,estCtrl)
    summary(lT1)  
    plot(lT1)
    
    # Check that the LL is higher than the target:
    lT1$Estimated$value > ll.target 
    
  }
  
  # Filter out TV & estimate GARCH component
  if(TRUE) {
    e <- e_orig/sqrt(lT1@g)
    
    lG1 <- MTV1@garch                # For first iteration
    lG1 <- MTV1$results[[1]]$garch   # For subsequent iterations - update the index!
    summary(lG1)
    ll.target <- getTargetValue(e,garchObj = lG1)
    lG1$pars <- lG1$Estimated$pars
    
    lG1$pars["omega",1] <- 0.10
    lG1$pars["alpha",1] <- 0.05
    lG1$pars["beta",1] <- 0.80
    lG1$pars["gamma",1] <- 0.02
    
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
  MTV1$Estimated[[1]] <- MTV1$results[[1]]
  
  
  # Save the estimated model
  saveRDS(MTV1,"Output/ANZ_mtvgjr_manual.RDS")
}

##==============================================================##
##                            THE END
##==============================================================##


summary(MTV1$results[[1]]$tv)
summary(MTV1$results[[2]]$tv)

