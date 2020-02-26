
#### ====================== Initialisation ======================== ####
rm(list=ls())
gc(T)

library(graphics)
library(foreach)
#library(doParallel)
## Note the RevoUtilsMath library parallelises the matrix calculations, so we don't need to use doParallel anymore
library(RevoUtilsMath)
library(moments)
library(stats)
library(Matrix)
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
  #e_cba <- mydata$e_cba
  e_nab <- mydata$e_nab
  #e_wbc <- mydata$e_wbc
  rm(mydata)
}
##====================== Data Load =============================##


#### ====================== Plot the Data ========================= ####
DoThis <- TRUE
if (DoThis) {
  ptitle <- "NAB de-meaned Returns"
  plot(dates,e_nab,"l",main = ptitle)
  plot(e_nab,type="l",main = ptitle)
  abline(v=seq(1,NROW(e_nab),by = 100),lwd=0.25,col="lightgrey")
  lines(e_nab,type="l")
  
  # Now smooth the squared returns:
  ptitle <- "NAB Squared Returns"
  plot(sqrt(e_nab^2),type="l",main = ptitle)
  abline(v=seq(1,NROW(e_nab),by = 100),lwd=0.25,col="lightgrey")
  lines(sqrt(e_nab^2),type="l")
  e2s <- ksmooth(seq(1:NROW(e_nab)),e_nab^2,kernel = "normal",bandwidth = 1100)    
  ptitle <- "NAB Squared Returns - Smoothed"
  plot(e2s,type="l",main = ptitle)
  abline(v=seq(1,NROW(e_nab),by = 100),lwd=0.25,col="lightgrey")
}
##====================== Plot the Data =========================##


####=================== Remove 'g' & estimate Garch ======================####
  
  if(FALSE) {
    e_std <- e_nab/sqrt(e2s$y)
  
  # Calculate the BIC
  lm1 <- lm(e_std ~ 1)
  #summary(lm1)
  bic1 <- BIC(lm1)
  
  ## BIC values for lm(e_std ~ 1), based on different Bandwidth settings: (Optimum = approx. 1100)
  #200 = 19686
  #300 = 19686
  #400 = 19626
  #500 = 19579
  #600 = 19527
  #700 = 19480
  #800 = 19443
  #900 = 19416
  #1000 = 19401
    #1100 = 19396 **
  #1200 = 19398
  #1300 = 19407
  #1500 = 19433
  #2000 = 19484
  
  ptitle <- "NAB Standardised Returns"
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
##   TV - AUS_4Banks\NAB ---- ORDER 0 
##==============================================================##

##  Method: Search across the data series, by changing the subset in the first line below
##          Log the results for each test below the code.

  e <- e_nab[5400:7027]
  # Create a TV object with just delta0 & st:
  ST <- 1:length(e)/length(e)
  TV <- newTV(st = ST,del0 = 1)
  # Estimate the univar series:
  TV <- EstimateTV(e,TV)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0

  if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_2018.RDS")
  fName <- "TestStatDist_NAB.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1100)
  
  TEST$FailCount_TR2
  TEST$FailCount_Robust
  
  hist(TEST$Stat_TR2,30)
  hist(TEST$Stat_Robust,30)

#### ========  SPECIFICATION - ORDER 0  - Results  ======== ####  
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 1:3400            ##
  
  ## Result: Pvals=[TR2: 15.4%, Robust: 3.45%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -5488.862
  ## delta0: 1.478184
  ## Pars: NA
  #
  ## TR2: Fail to Reject the null: H0=Model is sufficient
  ## Robust: Reject the null: H0=Model is suspect
  ##==============================================================##
  
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 1:2000            ##
  
  ## Result: Pvals=[TR2: 9.73%, Robust: 8.64%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -3118.776
  ## delta0: 1.323272
  ## Pars: NA
  #
  ## TR2: Fail to Reject the null: H0=Model is sufficient
  ## Robust: Fail to Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 2001:3600            ##
  
  ## Result: Pvals=[TR2: 43.6%, Robust: 1.82%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -2625.407
  ## delta0: 1.559721
  ## Pars: NA
  #
  ## TR2: Fail to Reject the null: H0=Model is sufficient
  ## Robust: Reject the null: H0=Model is suspect
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 1:4000            ##
  
  ## Result: Pvals=[TR2: 9.09%, Robust: 1.45%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -6317.862
  ## delta0: 1.378217
  ## Pars: NA
  #
  ## TR2: Fail to Reject the null @ 5%: H0=Model is sufficient
  ## Robust: Reject the null: H0=Model is suspect
  ##   We will estimate a 1-Order model over this range, since
  ##   the rejection is stronger than for 1:3400
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 3200:4400         ##
  
  ## Result: Pvals=[TR2: 0.09%, Robust: 1.0%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -6317.862
  ## delta0: 1.378217
  ## Pars: NA
  #
  ## TR2: Reject the null: H0=Model is suspect
  ## Robust: Reject the null: H0=Model is suspect
  ##   We will estimate a 1-Order model over this range
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 4200:5600         ##
  
  ## Result: Pvals=[TR2: 0.64%, Robust: 0.46%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -2995.351
  ## delta0: 4.215567
  ## Pars: NA
  #
  ## TR2: Reject the null: H0=Model is suspect
  ## Robust: Reject the null: H0=Model is suspect
  ##   We will estimate a 1-Order model over this range
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 0 ---- obs 5400:7027         ##
  
  ## Result: Pvals=[TR2: 12.5%, Robust: 21.5%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -2467.551
  ## delta0: 1.213779
  ## Pars: NA
  #
  ## TR2: Fail to Reject the null: H0=Model is sufficient
  ## Robust: Fail to Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  
  
  
  
#### ---  H3,H2,H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H3 
##==============================================================##

  e <- e_nab[1:2500]
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

  if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_2018.RDS")
  fName <- "TestStatDist_NAB.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,testorder = TESTorder$H03, saveas = fName,numloops = 1100)

 
  ##  Could not calculate RefTest stats - error in Solve()
  
##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H3 

## Result: Pvals=[TR2: 25.8%, Robust: 26%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -3130.408
## delta0: 0.9148737
## Pars: NA
## Linpars: -0.09732136  1.62535678
#
## Fail to Reject the NULL.  No evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H2 
##==============================================================##

  e <- e_nab[1:4000]
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
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 1100)
  
  ##  Could not calculate RefTest stats - error in Solve()
  
##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H2 

## Result: Pvals=[TR2: 32.18%, Robust: 36.55%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -3133.77
## delta0: 0.7064782
## Pars: NA
## Linpars: 1.397693
#
## Cannot Reject the NULL.  No evidence for quadratic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 0 /H1 
##==============================================================##

  e <- e_nab[1:4000]
  # Create a TV object with just delta0 & st:
  ST <- 1:length(e)/length(e)
  TV <- newTV(del0 = 1,st = ST)
  # Estimate the univar series:
  TV <- EstimateTV(e,TV)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$linpars <- c(NA)
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_2018.RDS")
  fName <- "TestStatDist_NAB.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,testorder = TESTorder$H01, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 0 H1 

## Result: Pvals=[TR2: 0.818%, Robust: 0.182%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -6317.862
## delta0: 1.378217
## Pars: NA
## Linpars: NA
#
## Reject the NULL.  Strongest evidence for linear transition
## Now, let's estimate this transition...
##==============================================================##


  
  
  
#### ========  SPECIFICATION - ORDER 1  ======== ####
##==============================================================##
##   TV - AUS_4Banks\NAB ---- ORDER 1 
##==============================================================##
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
  
  e <- e_nab[4200:5600]
  # Create a TV object with just delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,3,0.5)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$single,speedopt=TRspeedopt$eta)
  startpars <- c(1,3,0.2,0.7)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$double,speedopt=TRspeedopt$eta)
  
  # Estimate the univar series:
  parScale <- c(2,2,1,0.5) # Single
  parScale <- c(1,1,1,0.5,0.5) # Double
  nDeps <- rep(1e-5,length(parScale))
  nDeps <- c(1e-5,1e-5,1e-5,1e-7,1e-7)
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7)
  TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$pars
  plot(TV$Estimated$condvars,type = 'l')

  if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_2018.RDS")
  fName <- "TestStatDist_NAB.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1200)
  
#### ========  SPECIFICATION - ORDER 1  - Results  ======== ####
##==============================================================##
##  TV - AUS_4Banks\NAB ---- ORDER 1 --- Obs: 1:4000

## Result: Pvals=[TR2: 82.3%, Robust: 74.9%], using
## TestStat_ProbDist_NAB Logliklihood Value:  -6143.415
## delta0: 2.263151
## Pars: 
##             [,1]
##  delta -1.3706683
##  speed  4.9439988
##  loc1   0.3432577
##  loc2   0.7063913
#
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
##==============================================================##

  
  #==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 1 --- Obs: 3200:4400
  
  ## Result: Pvals=[TR2: 80.7%, Robust: 56.4%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -1954.609
  ## delta0: 0.7740263
  ## Pars: 
  ##             [,1]
  ##  delta 26.2780336
  ##  speed  2.8468008
  ##  loc1   0.9177653
  ##  loc2          NA
  #
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  #==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 1 --- Obs: 4200:5600
  
  ## Result: Pvals=[TR2: 43.5%, Robust: 25.8%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -2752.955
  ## delta0: 13.92481
  ## Pars: 
  ##             [,1]
  ##  delta -12.2182993
  ##  speed   2.8360982
  ##  loc1    0.2036849
  ##  loc2          NA
  #
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ## Note: These p-values are still relatively low...
  ##       
  ##==============================================================##
  
  
  
  #### ========  SPECIFICATION - ORDER 2    ======== ####
  
  ##==============================================================##
  ##   TV - AUS_4Banks\NAB ---- ORDER 2 
  ##==============================================================##
  
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
  
  
  e <- e_nab[1:4400]
  # Create a TV object with just delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,3,0.2,0.5,1,3,0.9)
   startpars <- c(1,3,0.2,0.5,1,3,0.9)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$double,TRshape$single),speedopt=TRspeedopt$eta)
  #startpars <- c(1,3,0.33,0.66)
  #TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$double,speedopt=TRspeedopt$eta)
  
  # Estimate the univar series:
  parScale <- c(1,1,1,0.7,0.7,1,1,0.7)
   parScale <- c(1,1,1,0.7,0.7,1,1,0.7)
  nDeps <- rep(1e-5,length(parScale))
   nDeps <- c(1e-5,1e-5,1e-5,1e-4,1e-4,1e-5,1e-5,1e-4)
  
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-9)
  TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$pars
  plot(TV$Estimated$condvars,type = 'l')
  TV <- calcStderr_TV(e,TV)
  TV$Estimated$hessian
  TV$Estimated$stderr
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_2018.RDS")
  fName <- "TestStatDist_NAB.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 2000)
  
  
  #### ========  SPECIFICATION - ORDER 2  - Results  ======== ####
  
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 2 ---- obs 1:4400            ##
  ##
  ## Result: Pvals=[TR2: 98.6%, Robust: 90.7%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -7048.411
  ## delta0: 2.263036
  ## Pars: 
  ##             [,1]       [,2]
  ##  delta -1.3777879 22.2423704
  ##  speed  5.1261674  4.2829338
  ##  loc1   0.3115018  0.9715571
  ##  loc2   0.6423641         NA
  ##
  ##  StdErr: 0.117160440 0.120547651 0.229945182 0.010147392 0.009369571 
  ##          6.505067624 0.168930789 0.009459104
  ## 
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  
  
  #### ========  SPECIFICATION - ORDER 3    ======== ####
  
  ##==============================================================##
  ##   TV - AUS_4Banks\NAB ---- ORDER 3 
  ##==============================================================##
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
  
  e <- e_nab[1:5600]
  e <- e_nab[1:6149]
  e <- e_nab
  # Create a TV object with delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,3,0.2,0.5,1,3,0.7,1,3,0.9)  # [1:5600]
  startpars <- c(1,3,0.1,0.3,1,3,0.6,1,3,0.8)  # [1:7027]
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$double,TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)
  
  # Estimate the univar series:
  startpars <- c(1,3,0.1,0.3, 1,3,0.6, 1,3,0.8)  # [1:7027]
  parScale <- c(3,3,3,0.05,0.05,5,3,0.05,5,3,0.05)
  nDeps <- rep(1e-5,length(parScale))
  
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
  
  TV <- calcParamStats_TV(TV)
  
  TV$Estimated$parsVector
  TV$Estimated$stderr
  TV$Estimated$tStat
  TV$Estimated$PValues
  
  saveRDS(TV,"Estimated_TVOrd3_NAB.RDS")
  
  TV <- readRDS("Estimated_TVOrd3_NAB.RDS")
  View(TV)
  
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_NAB_2018.RDS")
  fName <- "TestStatDist_NAB.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1600)
  
  
  #### ========  SPECIFICATION - ORDER 3  - Results  ======== ####
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 3 ---- obs 1:5600            ##
  ##
  ## Result: Pvals=[TR2: %, Robust: %], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -9253.566
  ## delta0: 2.263575* se: 0.116889052 
  ## Pars: 
  ##            [,1]       [,2]       [,3]
  ##  delta -1.3700961* 39.6228650* -38.7975651*  
  ##  speed  5.6168424*  4.4500735*   4.1803120*
  ##  loc1   0.2452802*  0.7618828*   0.7795491*
  ##  loc2   0.5045247*         NA          NA
  ## TR2 - 
  ## Robust - 
  ## Standard Errors:
  ##  [1] 0.120464362 0.230716440 0.007992835 0.007420259 9.591142235 0.088157739
  ##  [8] 0.003216139 9.588562996 0.096676734 0.003926781
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\NAB ---- ORDER 3 ---- obs 1:7027            ##
  ##
  ## Result: Pvals=[TR2: 78.5%, Robust: 54.9%], using
  ## TestStat_ProbDist_NAB Logliklihood Value:  -11438.12
  ## delta0: 2.263294
  ## Pars: 
  ##            [,1]       [,2]       [,3]
  ##  delta -1.3629961 23.1683258 -22.6972043
  ##  speed  6.0774728  4.7347456   4.2034498
  ##  loc1   0.1957340  0.6023369   0.6273989
  ##  loc2   0.4019763         NA          NA
  ##
  ##
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  #
  #
  # > TV$Estimated$parsVector
  # delta0     delta1     speed1     loc1.1     loc1.2     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
  # 2.263294  -1.362996   6.077473   0.195734   0.401976  23.168326   4.734746   0.602337 -22.697204   4.203450   0.627399 
  # > TV$Estimated$stderr
  # delta0   delta1   speed1   loc1.1   loc1.2   delta2   speed2   loc2.1   delta3   speed3   loc3.1 
  # 0.116507 0.120179 0.231195 0.006381 0.005945 7.742008 0.158801 0.003530 7.736216 0.093737 0.005604 
  # > TV$Estimated$tStat
  # delta0     delta1     speed1     loc1.1     loc1.2     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
  # 19.426250 -11.341382  26.287216  30.674502  67.615812   2.992547  29.815593 170.633711  -2.933890  44.843018 111.955567 
  # > TV$Estimated$PValues
  # delta0*  delta1*  speed1*  loc1.1*  loc1.2*  delta2*  speed2*  loc2.1*  delta3*  speed3*  loc3.1* 
  #   0.000000 0.000000 0.000000 0.000000 0.000000 0.002637 0.000000 0.000000 0.003190 0.000000 0.000000 
  # 
  
#### ====================  Conclusion  ======================== ####
   
##
## Recommend the Ord3 model
## 
##==============================================================##


##==============================================================##
##                            THE END
##==============================================================##


