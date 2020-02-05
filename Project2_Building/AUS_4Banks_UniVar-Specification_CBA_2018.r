
#### ====================== Initialisation ======================== ####
rm(list=ls())
gc(T)

library(graphics)
library(foreach)
#library(doParallel)
## Note the RevoUtilsMath library parallelises the matrix calculations, so we don't need to use doParallel anymore
library(RevoUtilsMath)
library(moments)
library(Matrix)
library(stats)
library(stats4)


## References: List of other XXX.R files used:
## 1. AUS_4Banks-GenRefData.R

#setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC & Laptop - GOOGLE DRIVE
#setwd("E:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Glen's Laptop - GOOGLE DRIVE
setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Glen's home PC - GOOGLE DRIVE

functionsPath <- file.path(dirname(getwd()),"Functions")
functionsFile <- file.path(functionsPath,"functions_tvgjr_v7.r")
source(functionsFile)

# ## Setup the parallel backend environment ##
# numcores <- 3
# Sys.setenv("MC_CORES" = numcores)
# cl <- makeCluster(numcores)
# registerDoParallel(cl, cores = numcores)

# # Remove parallel backend - when you need to ##
# unregisterDoParallel(cl)
# rm(cl)

##====================== Initialisation ========================##


#### ====================== Data Setup ============================ ####
DoThis <- FALSE
if (DoThis){
  mydata <- read.csv("Bank-Returns-2018DEC07.csv",header=TRUE)
  dates <- as.Date(mydata$date, format = "%d/%m/%Y")
  
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
  saveRDS(mydata,"AUS_4Banks-ReturnsData-2018.RDS")
  saveRDS(dates,"AUS_4Banks-Dates-2018.RDS")
}
##====================== Data Setup ============================##

#### ====================== Data Load ============================= ####
if (TRUE){
  # Read AUS_4Banks data from saved file
  dates <- readRDS("AUS_4Banks-Dates-2018.RDS")
  mydata <- readRDS("AUS_4Banks-ReturnsData-2018.RDS")
  #e_anz <- mydata$e_anz
  e_cba <- mydata$e_cba
  #e_nab <- mydata$e_nab
  #e_wbc <- mydata$e_wbc
  rm(mydata)
}
##====================== Data Load =============================##


#### ====================== Plot the Data ========================= ####
DoThis <- TRUE
if (DoThis) {
  ptitle <- "CBA de-meaned Returns"
  plot(dates,e_cba,"l",main = ptitle)
  plot(e_cba,type="l",main = ptitle)
  abline(v=seq(1,NROW(e_cba),by = 100),lwd=0.25,col="lightgrey")
  lines(e_cba,type="l")
  
  # Now smooth the squared returns:
  ptitle <- "CBA Squared Returns"
  plot(sqrt(e_cba^2),type="l",main = ptitle)
  abline(v=seq(1,NROW(e_cba),by = 100),lwd=0.25,col="lightgrey")
  lines(sqrt(e_cba^2),type="l")
  e2s <- ksmooth(seq(1:NROW(e_cba)),e_cba^2,kernel = "normal",bandwidth = 1100)    
  ptitle <- "CBA Squared Returns - Smoothed"
  plot(e2s,type="l",main = ptitle)
  abline(v=seq(1,NROW(e_cba),by = 100),lwd=0.2,col="lightgrey")
  
}
##====================== Plot the Data =========================##

####=================== Remove 'g' & estimate Garch ======================####

if (TRUE) {
  
  e_std <- e_cba/sqrt(e2s$y)
  
  # Calculate the BIC
  lm1 <- lm(e_std ~ 1)
  #summary(lm1)
  bic1 <- BIC(lm1)
  
  ## BIC values for lm(e_std ~ 1), based on different Bandwidth settings: (Optimum = approx. 1100)
  
  #800 = 19525
  #900 = 19501
  #1000 = 19488  
  #1100 = 19486 **
  #1200 = 19491
  #1300 = 19504
  #1500 = 19537

  ptitle <- "CBA Standardised Returns"
  plot(e_std,type="l",main = ptitle)
  G1 <- newGARCH(GARCHtype$general,c(0.05,0.05,0.85))
  G1 <- EstimateGARCH(e_std,G1)
  G1$Estimated$pars
  ##  Use these parameters to generate a reference data set for use in Bootstrapping
  ##  a distribution when testing.  See AUS_4Banks-GenRefData.R
}

##=================== Remove 'g' & estimate Garch ======================##
  
##=========================================================================================================================##
##   General Specification Method:
##
##   1. Kernel-Smooth the demeaned returns and identify the Bandwith that optimises BIC
##   2. Divide this smoothed series out of the returns data & then estimate the remaining GARCH
##   2.1  Use these Garch Parameters to create a reference-data set that will be used for bootstrapping 
##        a reference distribution during testing
##   3. Iteratively test sub-sections of the sample for a transition (H0 Test)
##   3.1  Experience has shown that the test can fail to identify transitions when there are too many in a sample
##        Further research is required to fully understand why.  Best guess is that the sum of all transitions
##        averages out to be close to a constant variance (i.e. no transition)
##        Therefore we recommend testing subsections.  Use the smoothed square returns plot as a guide
##        for the sub-section intervals to test.
##   3.2  Important!  Each sub-section should be a minimum of 1000 observations
##   4. Once evidence of a transition is found run the 3 aternate hypothesis tests to identify the *most likely* shape
##   4.1  Note: Visual observation of the smoothed square returns can give clues to the shape
##   5. Now re-estimate the sub-sample with the identified Order 1 model, and repeat the H0 Test
##   5.1  Note: Failure to reject the NULL here indicates there are no more transitions in this sub-section
##   6. Repeat steps 3 - 5 for each subsequent sub-section
##   6.1  Note: It is ok to allow small overlaps in your subsections when testing
##        E.g. Subsection#2 = [2000:3200], Subsection#3 = [2900:4500]
##   7. Complete the specification by adding the individual sub-sections together
##      i.e.  If we have identified 3 individual transitions in our sample, then estimate an Ord 3 model
##  8. Consider if any of the transitions could/should be parameterised as squares or cubes?
##  9.1  This can be tested by estimating the alternat model and comparing the LogLiklihood value & conditional variance plot
##  10. Estimate the final model and finally calculate & check the standard errors & parameter statistics
##  10.1  Note:  It is obviously ideal to define a model where all parameters are significant.
##               See the section on Optim-Tuning Tips & Tricks to optimise the parameterization
##               Very often the params in the first/most optimum model will not all be significant
##               It is often possible to refine this model (using optim-tuning) to get a model
##               with a very similar Logliklihood (usually slightly lower), but with all params being significant
##        Essentially, we can trade-off our 'best numerical fit' (logliklihood value) against parameter-significance.
##
##=========================================================================================================================##


#### ========  SPECIFICATION - ORDER 0  ======== ####
##==============================================================##
##   TV - AUS_4Banks\CBA ---- ORDER 0 
##==============================================================##

##  Method: Search across the data series, by changing the subset in the first line below
##          Log the results for each test below the code.

  e <- e_cba[6400:7027]
  # Create a TV object with just delta0 & st:
  ST <- 1:length(e)/length(e)
  TV <- newTV(st = ST,del0 = 1)
  # Estimate the univar series:
  TV <- EstimateTV(e,TV)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0

  if (!exists("refData")) refData <- readRDS("RefData_withGarch_CBA_2018.RDS")
  fName <- "TestStatDist_CBA.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 3000)
  
  TEST$FailCount_TR2
  TEST$FailCount_Robust
  
  hist(TEST$Stat_TR2,30)
  hist(TEST$Stat_Robust,30)

#### ========  SPECIFICATION - ORDER 0  - Results  ======== ####  
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 1:3300            ##
  ##
  ## Result: Pvals=[TR2: 27.2%, Robust: 7.09%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -5028.909
  ## delta0: 1.233836
  ## Pars: NA
  #
  ## Fail to Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 400:2000          ##
  ##
  ## Result: Pvals=[TR2: 54%, Robust: 49.7%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -2444.777
  ## delta0: 1.240557
  ## Pars: NA
  #
  ## Fail to Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 1:4000          ##
  ##
  ## Result: Pvals=[TR2: 13.1%, Robust: 6.55%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -5928.433
  ## delta0: 1.134587
  ## Pars: NA
  #
  ## Fail to Reject the null: H0=Model is sufficient
  ## Note: These p-values are quite low.  Remember that the reference distribution
  ##       is just a boot-strapped approximation.  There is a possibility of
  ##       a transition here.
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 2500:3500          ##
  ##
  ## Result: Pvals=[TR2: 14.3%, Robust: 3.55%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -1460.564
  ## delta0: 1.084701
  ## Pars: NA
  #
  ##    TR2: Fail to Reject the null: H0=Model is sufficient
  ## Robust: Reject the null: H0=Model is suspect
  ##
  ## Now let's check an Order1 model & compare LogLik values...
  ## See below: Order1 model produces LogLik = -1401.763
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 1:2500          ##
  ##
  ## Result: Pvals=[TR2: 57.9%, Robust: 62.1%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -3829.767
  ## delta0: 1.253321
  ## Pars: NA
  #
  ## Fail to Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 3500:4400         ##
  ##
  ## Result: Pvals=[TR2: 0.273%, Robust: 0.909%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -1708.805
  ## delta0: 2.600903
  ## Pars: NA
  #
  ## STRONGLY Reject the null: H0=Model is suspect
  ## Now let's check an Order1 model & compare LogLik values...
  ## See below: Order1 model produces LogLik = -1437.55
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 4200:5400         ##
  ##
  ## Result: Pvals=[TR2: 1.45%, Robust: 1.09%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -2464.63
  ## delta0: 3.550774
  ## Pars: NA
  #
  ## STRONGLY Reject the null: H0=Model is suspect
  ## Now let's check an Order1 model & compare LogLik values...
  ## See below: Order1 model produces LogLik = -2253.822
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 5400:6400         ##
  ##
  ## Result: Pvals=[TR2: 10.7%, Robust: 4.7%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -1501.471
  ## delta0: 1.177111
  ## numloops for Simulatied Dist = 3000
  ## Pars: NA
  #
  ##    TR2: Fail to Reject the null: H0=Model is sufficient
  ## Robust: Reject the null @ 5%: H0=Model is suspect
  ##
  ## Now let's check an Order1 model & compare LogLik values...
  ## See below: Order1 model produces LogLik = -1437.355
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 5700:6700         ##
  ##
  ## Result: Pvals=[TR2: 18.9%, Robust: 4.83%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -1516.149
  ## delta0: 1.210346
  ## numloops for Simulatied Dist = 3000
  ## Pars: NA
  #
  ##    TR2: Fail to Reject the null: H0=Model is sufficient
  ## Robust: Reject the null @ 5%: H0=Model is suspect
  ##
  ## Validates the above result for [5400:6400]
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 0 ---- obs 6400:7027         ##
  ##
  ## Result: Pvals=[TR2: 97.8%, Robust: 97.6%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -877.6372
  ## delta0: 0.9559591
  ## numloops for Simulatied Dist = 3000
  ## Pars: NA
  #
  ##    TR2: Fail to Reject the null: H0=Model is sufficient
  ## Robust: Fail to Reject the null: H0=Model is sufficient
  ##
  ##==============================================================##
  
#### ---  H3,H2,H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\CBA ---- ORDER 0 /H3 
##==============================================================##

  e <- e_cba[1:4000]
  # Create a TV object with just delta0 & st:
  ST <- (1:length(e))/length(e)
  TV <- newTV(st = ST,del0 = 1)
  # Estimate the univar series:
  TV <- addH03pars(TV)
  TV <- EstimateTV(e,TV)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$linpars
  plot(TV$Estimated$condvars,type = 'l')
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_CBA_2018.RDS")
  fName <- "TestStatDist_CBA.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H03)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H03)
  TEST <- GenTestStatDist(e,TV,refData,refTest,testorder = TESTorder$H03, saveas = fName,numloops = 1100)

 
  ##  Could not calculate RefTest stats - error in Solve()
  
##==============================================================##
##  TV - AUS_4Banks\CBA ---- ORDER 0 H3 -- obs 1:4000

## Result: Pvals=[TR2: %, Robust: %], using
## TestStat_ProbDist_CBA Logliklihood Value:  -6731.278
## delta0: 2.542194
## Pars: NA
## Linpars: -0.4627832 -1.5030424
#
## Fail to Reject the NULL.  No evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\CBA ---- ORDER 0 /H2 
##==============================================================##

  e <- e_cba[1:4000]
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
  plot(TV$Estimated$condvars,type = 'l')
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_CBA_2018.RDS")
  fName <- "TestStatDist_CBA.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H02)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H02)
  TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 1100)
  
  ##  Could not calculate RefTest stats - error in Solve()
  
##==============================================================##
##  TV - AUS_4Banks\CBA ---- ORDER 0 H2 

## Result: Pvals=[TR2: 32.18%, Robust: 36.55%], using
## TestStat_ProbDist_CBA Logliklihood Value:  -3133.77
## delta0: 0.7064782
## Pars: NA
## Linpars: 1.397693
#
## Cannot Reject the NULL.  No evidence for quadratic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\CBA ---- ORDER 0 /H1 
##==============================================================##

  e <- e_cba[1:4000]
  # Create a TV object with just delta0 & st:
  ST <- 1:length(e)/length(e)
  TV <- newTV(del0 = 1,st = ST)
  # Estimate the univar series:
  TV <- EstimateTV(e,TV)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  plot(TV$Estimated$condvars,type = 'l')
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_CBA_2018.RDS")
  fName <- "TestStatDist_CBA.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H01)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H01)
  TEST <- GenTestStatDist(e,TV,refData,refTest,testorder = TESTorder$H01, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\CBA ---- ORDER 0 H1 

## Result: Pvals=[TR2: 0.818%, Robust: 0.182%], using
## TestStat_ProbDist_CBA Logliklihood Value:  -6317.862
## delta0: 1.378217
## Pars: NA
## Linpars: NA
#
## Reject the NULL.  Strongest evidence for linear transition
## Now, let's estimate this transition...
##==============================================================##

  
  
#### ========  SPECIFICATION - ORDER 1  ======== ####
##==============================================================##
##   TV - AUS_4Banks\CBA ---- ORDER 1 
##==============================================================##

  
##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
##          Log the results for each test below the code.
  
  e <- e_cba[5400:6400]
  # Create a TV object with just delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,3,0.5)
   startpars <- c(1,1,0.5)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$single,speedopt=TRspeedopt$eta)

  # Estimate the univar series:
  parScale <- c(1,1,1,0.5)
   parScale <- c(1,1,1,0.5)
  nDeps <- rep(1e-5,length(parScale))
   nDeps <- c(1e-5,1e-5,1e-5,1e-7)
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7)
  TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$pars  
  plot(TV$Estimated$condvars,type = 'l')
  TV <- calcStderr_TV(e,TV)
  TV$Estimated$hessian
  TV$Estimated$stderr

  if (!exists("refData")) refData <- readRDS("RefData_withGarch_CBA_2018.RDS")
  fName <- "TestStatDist_CBA.RDS"
  # Calculate the Reference test statistics:
  
  ##  Note: The optimiser can be quite sensitive to starting values and optim parameters
  ##        The parameter results returned may be a local maximum, or may have gotten
  ##        stuck on a boundary.
  ##        It is worth checking that the estimated model looks 'reasonable'
  ##        If the Tests below cannot be calculated, this is often a result of a
  ##        poorly estimated model.
  
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 3000)
  
  # Confirm we managed to execute more than 1000 samples for our distribution:
  TEST$FailCount_TR2
  TEST$FailCount_Robust
 

#### ========  SPECIFICATION - ORDER 1  - Results  ========
##==============================================================##
##  TV - AUS_4Banks\CBA ---- ORDER 1 ---- obs: 2500:3500

## Result: Pvals=[TR2: 64.3%, Robust: 68.2%], using
## TestStat_ProbDist_CBA Logliklihood Value:  -1401.763
## delta0: 1.649649
## Pars: 
##           [,1]
##  delta -1.0431751
##  speed  4.5619153
##  loc1   0.4579069
##  loc2         NA
#
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
##==============================================================##

  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 1 ---- obs: 3500:4400
  
  ## Result: Pvals=[TR2: 93%, Robust: 88.4%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -1437.55
  ## delta0: 0.6649733
  ## Pars: 
  ##           [,1]
  ##  delta 8.0653366
  ##  speed 3.1313321
  ##  loc1  0.7613485
  ##  loc2         NA
  #
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 1 ---- obs: 4200:5400
  
  ## Result: Pvals=[TR2: 85%, Robust: 62.7%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -2253.822
  ## delta0: 10.39315
  ## Pars: 
  ##           [,1]
  ##  delta -9.1952307
  ##  speed  2.5252518
  ##  loc1   0.2526861
  ##  loc2         NA
  #
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 1 ---- obs: 5400:6400
  
  ## Result: Pvals=[TR2: 34.6%, Robust: 27.6%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -1437.355
  ## delta0: 0.7297468
  ## Pars: 
  ##           [,1]
  ##  delta 1.3783038
  ##  speed 4.4973561
  ##  loc1  0.6771262
  ##  loc2         NA
  #
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##  We have identified 4 transitions.
  ##  The middle 2 could be approximated by a double
  ##  The final transition seems to be barely significant, based on
  ##  the fact that the 1st order model LogLik value is only 
  ##  slightly better than the null no-transition model.
  
  
    
#### ========  SPECIFICATION - ORDER 2    ======== ####
  
  ##==============================================================##
  ##   TV - AUS_4Banks\CBA ---- ORDER 2 
  ##==============================================================##
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
  
  e <- e_cba[1:4500]
  # Create a TV object with just delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,3,0.5,3,3,0.7)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)
  #startpars <- c(1,3,0.33,0.66)
  #TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$double,speedopt=TRspeedopt$eta)
  
  # Estimate the univar series:
  parScale <- c(1,1,1,0.3,2,1,0.3)
  nDeps <- rep(1e-5,length(parScale))
  #nDeps <- c(1e-5,1e-5,1e-5,1e-9,1e-5,1e-5,1e-9)
  
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7)
  TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$pars
  plot(TV$Estimated$condvars,type = 'l')
  TV <- calcStderr_TV(e,TV)
  TV$Estimated$hessian
  TV$Estimated$stderr
  
  loc_idx <- round(TV$Estimated$pars["loc1",1] * length(e) + 4400)
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_CBA_2018.RDS")
  fName <- "TestStatDist_CBA.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1600)
  

#### ========  SPECIFICATION - ORDER 2  - Results  ======== ####

##==============================================================##
##  TV - AUS_4Banks\CBA ---- ORDER 2 ---- obs 1:4500            ##
##
## Result: Pvals=[TR2: 63.5%, Robust: 55.2%], using
## TestStat_ProbDist_CBA Logliklihood Value:  -6932.005
## delta0: 1.310307
## Pars: 
##             [,1]       [,2]
##  delta -0.6813599 9.2164982
##  speed  6.3857010 4.6330189
##  loc1   0.6590038 0.9325437
##  loc2          NA         NA
## 
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
##==============================================================##


  
#### ========  SPECIFICATION - ORDER 3    ======== ####
  
  ##==============================================================##
  ##   TV - AUS_4Banks\CBA ---- ORDER 3 
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
  
 
#### ========  SPECIFICATION - ORDER 3  - Results  ======== ####
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 3 ---- obs 1:5400            ##
  ##
  ## Result: Pvals=[TR2: 75.1%, Robust: 49.3%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -8426.928
  ## delta0: 1.310056
  ## Pars: 
  ##            [,1]       [,2]       [,3]
  ##  delta -0.6805553 14.1947237 -13.594348
  ##  speed  6.5706915  4.6001632   4.067647
  ##  loc1   0.5491619  0.7825051   0.824174
  ##  loc2          NA         NA         NA
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 3 ---- obs 1:7027            ##
  ##
  ## Result: Pvals=[TR2: %, Robust: %], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -10817.32 - Can't calc StdErrs
  ## delta0: 1.311469
  ## Pars: 
  ##            [,1]       [,2]       [,3]
  ##  delta -0.6131891 8.508880 -8.0821080
  ##  speed  6.9999999 6.152717  4.4694729
  ##  loc1   0.4215623 0.595472  0.6474193
  ##  loc2          NA         NA         NA
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
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
  
  
  
  
#### ========  SPECIFICATION - ORDER 4    ======== ####
  
  ##==============================================================##
  ##   TV - AUS_4Banks\CBA ---- ORDER 4 
  ##==============================================================##
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
  
  e <- e_cba
  # Create a TV object with delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,3,0.3,3,3,0.5,1,3,0.7,1,3,0.9)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single,TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)
  
  # Estimate the univar series:
  parScale <- c(1,3,1,0.3,3,1,0.3,3,1,0.3,1,1,0.3)
  # Estimate the univar series:
  nDeps <- c(1e-5,1e-5,1e-7,1e-9,1e-5,1e-7,1e-9,1e-5,1e-7,1e-9,1e-5,1e-7,1e-9)  #LL = -10778.22
  nDeps <- c(1e-5,1e-5,1e-5,1e-7,1e-5,1e-5,1e-7,1e-5,1e-5,1e-7,1e-5,1e-5,1e-7)  #LL = -10777.77
  
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7,maxit=1000)
  # To get a smoother model where we can calc StdErr's
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-6,maxit=1000)
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
  
  saveRDS(TV,"Estimated_TVOrd4_CBA.RDS")
  
  TV <- readRDS("Estimated_TVOrd4_CBA.RDS")
  View(TV)
  
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_CBA_2018.RDS")
  fName <- "TestStatDist_CBA.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 5000)
  #
  # LMTR2 ref stat could not be calculated, so I will set it = Robust
  refTest$LMTR2 <- refTest$LMRobust
  # We will ignore the TR2 dist result (code bugs force us to run both atm)...
  
  
  
  
#### ========  SPECIFICATION - ORDER 4  - Results  ======== ####
  
  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 4 ---- obs 1:7027            ##
  ##
  ## Result: Pvals=[TR2: ---%, Robust: 92.5%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -10777.79 - Can't Calc StdErrs
  ## delta0: 3.345083e-06
  ## Pars: 
  ##            [,1]       [,2]       [,3]       [,4]
  ##  delta -1.1398010 11.1963760 -11.3993972 4.1994036
  ##  speed  5.7911382  4.9866985   3.9960316 0.6130525
  ##  loc1   0.4211042  0.5973416   0.6393683 0.6541467
  ##  loc2          NA         NA         NA      NA
  ##
  ## TR2 - Could not calculate Ref Test Stat
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##
  ## Note: This model is only marginally better than the 3rd-Order
  ##       It may not be worth the extra paramter...
  ##==============================================================##
  

  ##==============================================================##
  ##  TV - AUS_4Banks\CBA ---- ORDER 4 ---- obs 1:7027            ##
  ##
  ## Result: Pvals=[TR2: ---%, Robust: 92.5%], using
  ## TestStat_ProbDist_CBA Logliklihood Value:  -10777.83 - Can Calc StdErrs!
  ## delta0: 0.000801408
  ## Pars: 
  ##            [,1]       [,2]       [,3]       [,4]
  ##  delta  -1.1800028 12.1069323 -12.4447117 4.6652081
  ##  speed  6.2199613  4.8974647   3.9391073 0.6442248
  ##  loc1   0.4210966  0.5976255   0.6366318 0.7198601
  ##  loc2          NA         NA         NA      NA
  ##
  ## TR2 - Could not calculate Ref Test Stat
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##
  ## Note: This model is only marginally better than the 3rd-Order
  ##       It may not be worth the extra paramter...
  ##
  ##
  # TV$Estimated$parsVector
  # delta0        delta1        speed1        loc1.1        delta2        speed2        loc2.1        delta3        speed3        loc3.1        delta4 
  # 0.000801408  -1.180003000   6.219961000   0.421097000  12.106932000   4.897465000   0.597625000 -12.444712000   3.939107000   0.636632000   4.665208000 
  # speed4        loc4.1 
  # 0.644225000   0.719860000 
  # > TV$Estimated$stderr
  # delta0   delta1   speed1   loc1.1   delta2   speed2   loc2.1   delta3   speed3   loc3.1   delta4   speed4   loc4.1 
  # 1.494703 0.258540 0.688732 0.001702 6.558848 0.380505 0.006209 6.233154 0.194040 0.012694 4.868697 1.422702 0.645899 
  # > TV$Estimated$tStat
  # delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1     delta4     speed4     loc4.1 
  # 0.000536  -4.564102   9.031032 247.413043   1.845893  12.870961  96.251409  -1.996535  20.300490  50.152198   0.958205   0.452818   1.114509 
  # > TV$Estimated$PValues
  # delta0  delta1*  speed1*  loc1.1*   delta2  speed2*  loc2.1*  delta3*  speed3*  loc3.1*   delta4   speed4   loc4.1 
  # 0.949594 0.000005 0.000000 0.000000 0.061702 0.000000 0.000000 0.043619 0.000000 0.000000 0.321093 0.618159 0.251844 
  # 
  
  ##==============================================================##
  
####  =========  Conclusion  =========  ####
  
##  Recommend using the Order 3 model  

##==============================================================##
##                            THE END
##==============================================================##
  
  
  