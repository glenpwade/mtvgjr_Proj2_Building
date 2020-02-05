
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
#setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Glen's home PC - GOOGLE DRIVE
setwd("C:/MTVGJR_MGARCH/Project2_Building")        # Glen's Work laptop


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
  #e_nab <- mydata$e_nab
  e_wbc <- mydata$e_wbc
  rm(mydata)
}
##====================== Data Load =============================##


#### ====================== Plot the Data ========================= ####
DoThis <- TRUE
if (DoThis) {
  ptitle <- "WBC de-meaned Returns"
  plot(dates,e_wbc,"l",main = ptitle)
  plot(e_wbc,type="l",main = ptitle)
  abline(v=seq(1,NROW(e_wbc),by = 100),lwd=0.25,col="lightgrey")
  lines(e_wbc,type="l")
  
  # Now smooth the squared returns:
  ptitle <- "WBC Squared Returns"
  plot(sqrt(e_wbc^2),type="l",main = ptitle)
  abline(v=seq(1,NROW(e_wbc),by = 100),lwd=0.25,col="lightgrey")
  lines(sqrt(e_wbc^2),type="l")
  e2s <- ksmooth(seq(1:NROW(e_wbc)),e_wbc^2,kernel = "normal",bandwidth = 1000)    
  ptitle <- "WBC Squared Returns - Smoothed"
  plot(e2s,type="l",main = ptitle)
  abline(v=seq(1,NROW(e_wbc),by = 100),lwd=0.2,col="lightgrey")
}

  
##====================== Plot the Data =========================##

####=================== Remove 'g' & estimate Garch ======================####

if(FALSE) {
  e_std <- e_wbc/sqrt(e2s$y)
  
  # Calculate the BIC
  lm1 <- lm(e_std ~ 1)
  #summary(lm1)
  bic1 <- BIC(lm1)
  
  ## BIC values for lm(e_std ~ 1), based on different Bandwidth settings: (Optimum = approx. 1100)
  
  #500 = 19664
  #600 = 19640
  #700 = 19617
  #800 = 19599
  #900 = 19587
  #1000 = 19583  **
  #1100 = 19585
  #1200 = 19593
  #1300 = 19603
  #1500 = 19629
  #2000 = 19688
  
  
  ptitle <- "WBC Standardised Returns"
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
##   TV - AUS_4Banks\WBC ---- ORDER 0 
##==============================================================##

##  Method: Search across the data series, by changing the subset in the first line below
##          Log the results for each test below the code.

  e <- e_wbc[5800:7027]
  # Create a TV object with just delta0 & st:
  ST <- 1:length(e)/length(e)
  TV <- newTV(st = ST,del0 = 1)
  # Estimate the univar series:
  TV <- EstimateTV(e,TV)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0

  if (!exists("refData")) refData <- readRDS("RefData_withGarch_WBC_2018.RDS")
  fName <- "TestStatDist_WBC.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  getMKLthreads()
  setMKLthreads(4)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1200)
  
  hist(TEST$Stat_TR2,30)
  hist(TEST$Stat_Robust,30)

#### ========  SPECIFICATION - ORDER 0  - Results  ======== ####  
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 0 ---- obs 1:1000            ##
  
  ## Result: Pvals=[TR2: 83.4%, Robust: 26.1%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -1778.155
  ## delta0: 2.053174
  ## Pars: NA
  #
  ## Fail to Reject the null: H0=Model is sufficient
  ##==============================================================##
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 0 ---- obs 1:2000            ##
  
  ## Result: Pvals=[TR2: 71.5%, Robust: 65.2%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -3561.334
  ## delta0: 2.062416
  ## Pars: NA
  #
  ## Fail to Reject the null: H0=Model is sufficient
  ##==============================================================##
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 0 ---- obs 2000:3400         ##
  
  ## Result: Pvals=[TR2: 58.3%, Robust: 8.64%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -2177.948
  ## delta0: 1.311857
  ## Pars: NA
  #
  ## Fail to Reject the null: H0=Model is sufficient
  ##  BUT: The Robust test stat is low, so there may be something near here
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 0 ---- obs 1:3400            ##
  
  ## Result: Pvals=[TR2: 44.4%, Robust: 3.0%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -5778.491
  ## delta0: 1.753028
  ## Pars: NA
  #
  ## TR2:  Fail to Reject the null: H0=Model is sufficient
  ## Robust: Reject the null: H0=Model is suspect
  ##  We could estimate an Order 1 model here & compare LL values
  ##  or extend the sample range & check again...
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 0 ---- obs 1:3800            ##
  
  ## Result: Pvals=[TR2: 23.3%, Robust: 1.42%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -6355.075
  ## delta0: 1.660283
  ## Pars: NA
  #
  ## TR2:  Fail to Reject the null: H0=Model is sufficient
  ## Robust: Reject the null: H0=Model is suspect
  ##  We will estimate an Order 1 model here & compare LL values
  ##==============================================================##
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 0 ---- obs 3700:4700            ##
  
  ## Result: Pvals=[TR2: 2.2%, Robust: 0.17%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -2129.897
  ## delta0: 4.131459
  ## Pars: NA
  #
  ## TR2:  Reject the null: H0=Model is suspect
  ## Robust: Reject the null: H0=Model is suspect
  ##  We will estimate an Order 1 model here & compare LL values
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 0 ---- obs 4400:5900            ##
  
  ## Result: Pvals=[TR2: 0.33%, Robust: 1.92%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -2785.158
  ## delta0: 2.395825
  ## Pars: NA
  #
  ## TR2:  Reject the null: H0=Model is suspect
  ## Robust: Reject the null: H0=Model is suspect
  ##  We will estimate an Order 1 model here & compare LL values
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 0 ---- obs 5800:7027         ##
  
  ## Result: Pvals=[TR2: 21.6%, Robust: 14.8%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -1970.272
  ## delta0: 1.447382
  ## Pars: NA
  #
  ## TR2:  Fail to Reject the null: H0=Model is sufficient
  ## Robust: Fail to Reject the null: H0=Model is sufficient
  ##  
  ##==============================================================##
  
  
  
  
  

#### ---  H3,H2,H1 Shape Tests --- ####
##==============================================================##
##   TV - AUS_4Banks\WBC ---- ORDER 0 /H3 
##==============================================================##

  e <- e_wbc[1:2500]
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

  if (!exists("refData")) refData <- readRDS("RefData_withGarch_WBC_2018.RDS")
  fName <- "TestStatDist_WBC.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,testorder = TESTorder$H03, saveas = fName,numloops = 1100)

 
  ##  Could not calculate RefTest stats - error in Solve()
  
##==============================================================##
##  TV - AUS_4Banks\WBC ---- ORDER 0 H3 

## Result: Pvals=[TR2: 25.8%, Robust: 26%], using
## TestStat_ProbDist_WBC Logliklihood Value:  -3130.408
## delta0: 0.9148737
## Pars: NA
## Linpars: -0.09732136  1.62535678
#
## Fail to Reject the NULL.  No evidence for cubic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\WBC ---- ORDER 0 /H2 
##==============================================================##

  e <- e_wbc[1:4000]
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
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_WBC_1.RDS")
  fName <- "TestStatDist_WBC.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,testorder = TESTorder$H02, saveas = fName,numloops = 1100)
  
  ##  Could not calculate RefTest stats - error in Solve()
  
##==============================================================##
##  TV - AUS_4Banks\WBC ---- ORDER 0 H2 

## Result: Pvals=[TR2: 32.18%, Robust: 36.55%], using
## TestStat_ProbDist_WBC Logliklihood Value:  -3133.77
## delta0: 0.7064782
## Pars: NA
## Linpars: 1.397693
#
## Cannot Reject the NULL.  No evidence for quadratic transition
##==============================================================##


##==============================================================##
##   TV - AUS_4Banks\WBC ---- ORDER 0 /H1 
##==============================================================##

  e <- e_wbc[1:4000]
  # Create a TV object with just delta0 & st:
  ST <- 1:length(e)/length(e)
  TV <- newTV(del0 = 1,st = ST)
  # Estimate the univar series:
  TV <- EstimateTV(e,TV)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$linpars <- c(NA)
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_WBC_2018.RDS")
  fName <- "TestStatDist_WBC.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,testorder = TESTorder$H01, saveas = fName,numloops = 1100)

##==============================================================##
##  TV - AUS_4Banks\WBC ---- ORDER 0 H1 

## Result: Pvals=[TR2: 0.818%, Robust: 0.182%], using
## TestStat_ProbDist_WBC Logliklihood Value:  -6317.862
## delta0: 1.378217
## Pars: NA
## Linpars: NA
#
## Reject the NULL.  Strongest evidence for linear transition
## Now, let's estimate this transition...
##==============================================================##


  
  
  

  
#### ========  SPECIFICATION - ORDER 1  ======== ####
##==============================================================##
##   TV - AUS_4Banks\WBC ---- ORDER 1 
##==============================================================##
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
  
  e <- e_wbc[4400:5900]
  # Create a TV object with just delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,2,0.5)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$single,speedopt=TRspeedopt$eta)
  #startpars <- c(1,3,0.33,0.66)
  #TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$double,speedopt=TRspeedopt$eta)
  #startpars <- c(1,3,0.5)
  #TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$triple,speedopt=TRspeedopt$eta)
  
  # Estimate the univar series:
  parScale <- c(3,3,1,0.3)
  #parScale <- c(1,1,1,0.3,0.3)
  nDeps <- rep(1e-7,length(parScale))
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-9)
  TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$pars
  plot(TV$Estimated$condvars,type = 'l')
  TV <- calcStderr_TV(e,TV)
  TV$Estimated$stderr
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_WBC_2018.RDS")
  fName <- "TestStatDist_WBC.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1100)
  
  View(refTest$LMTR2)
  summary(refTest$LMTR2)
  
  
#### ==============  SPECIFICATION - Order 1 - Results  ========= ####
  
##==============================================================##
##  TV - AUS_4Banks\WBC ---- ORDER 1 ---- Obs: 1:3800

## Result: Pvals=[TR2: 89.9%, Robust: 90.2%], using
## TestStat_ProbDist_WBC Logliklihood Value:  -6259.467
## delta0: 2.000157
## Pars: 
##           [,1]
## delta -1.2562824
## speed  2.7344028
## loc1   0.7304873
## loc2         NA
##
##  This model is better than the Ord0, where LL= -6355.075
## TR2 - Cannot Reject the null: H0=Model is sufficient
## Robust - Cannot Reject the null: H0=Model is sufficient
##==============================================================##

  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 1 ---- Obs: 3700:4700
  
  ## Result: Pvals=[TR2: 22.2%, Robust: 12.2%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -1937.369
  ## delta0: 0.8447678
  ## Pars: 
  ##           [,1]
  ## delta  6.2104578
  ## speed  3.3400911
  ## loc1   0.4719298
  ## loc2         NA
  ##
  ##  This model is better than the Ord0, where LL= -2129.897
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 1 ---- Obs: 4400:5900
  
  ## Result: Pvals=[TR2: 8.04%, Robust: 78.3%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -2630.086
  ## delta0: 0.8447678
  ## Pars: 
  ##           [,1]
  ## delta  -10.60756
  ## speed   1.604203
  ## loc1    6.662361e-04
  ## loc2         NA
  ##
  ##  This model is better than the Ord0, where LL= -2785.158
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  


  
  
  
  #### ========  SPECIFICATION - ORDER 2    ======== ####
  
  ##==============================================================##
  ##   TV - AUS_4Banks\WBC ---- ORDER 2 
  ##==============================================================##
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
  
  e <- e_wbc[1:4700]
  # Create a TV object with just delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,3,0.5,1,3,0.7)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)
  #startpars <- c(1,3,0.33,0.66)
  #TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=TRshape$double,speedopt=TRspeedopt$eta)
  
  # Estimate the univar series:
  parScale <- c(1,1,1,0.5,2,1,0.5)
  nDeps <- rep(1e-5,length(parScale))
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7)
  TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
  # Estimation Outputs:
  TV$Estimated$value
  TV$Estimated$delta0
  TV$Estimated$pars
  plot(TV$Estimated$condvars,type = 'l')
  TV <- calcStderr_TV(e,TV)
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_WBC_2018.RDS")
  fName <- "TestStatDist_WBC.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1200)
  
  
  #### ========  SPECIFICATION - ORDER 2  - Results  ======== ####
  
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 2 ---- obs 1:4700            ##
  ##
  ## Result: Pvals=[TR2: 45.4%, Robust: 37.5%], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -8051.796
  ## delta0: 2.205168
  ## Pars: 
  ##             [,1]       [,2]
  ##  delta -1.2572205 6.3155955
  ##  speed  2.9522115 4.7995467
  ##  loc1   0.5910678 0.8869727
  ##  loc2          NA        NA
  ## 
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  
  
  
  #### ========  SPECIFICATION - ORDER 3    ======== ####
  
  ##==============================================================##
  ##   TV - AUS_4Banks\WBC ---- ORDER 3 
  ##==============================================================##
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
  
  e <- e_wbc
  # Create a TV object with just delta0, pars & st:
  ST <- 1:length(e)/length(e)
  startpars <- c(1,3,0.4,3,3,0.6,-3,3,0.8)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)

  parScale <- c(1,3,1,0.2,3,1,0.2,3,1,0.2)
  nDeps <- rep(1e-5,length(parScale))
  #nDeps <- c(1e-5,1e-5,1e-5,1e-9,1e-5,1e-5,1e-9,1e-5,1e-5,1e-9)  
  
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
  
  saveRDS(TV,"Estimated_TVOrd3_WBC.RDS")
  
  TV <- readRDS("Estimated_TVOrd3_WBC.RDS")
  View(TV)
  
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_WBC_2018.RDS")
  fName <- "TestStatDist_WBC.RDS"
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
  ##  TV - AUS_4Banks\WBC ---- ORDER 3 ---- obs 1:7027            ##
  ##
  ## Result: Pvals=[TR2: %, Robust: %], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -11831.75
  ## delta0: 2.307224
  ## Pars: 
  ##            [,1]       [,2]       [,3]
  ##  delta -3.4292611 15.4868180 -13.3671361
  ##  speed  1.7834176  4.6558082   5.5189624
  ##  loc1   0.5585688  0.6062854   0.6340601
  ##  loc2          NA         NA          NA
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##
  # > TV$Estimated$parsVector
  # delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
  # 2.307224  -3.429261   1.783418   0.558569  15.486818   4.655808   0.606285 -13.367136   5.518962   0.634060 
  # > TV$Estimated$stderr
  # delta0   delta1   speed1   loc1.1   delta2   speed2   loc2.1   delta3   speed3   loc3.1 
  # 0.166833 0.338872 0.173069 0.027516 3.061389 0.120580 0.003812 3.061534 0.371265 0.001900 
  # > TV$Estimated$tStat
  # delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1 
  # 13.829543 -10.119635  10.304665  20.299789   5.058755  38.611776 159.046432  -4.366156  14.865290 333.715789 
  # > TV$Estimated$PValues
  # delta0* delta1* speed1* loc1.1* delta2* speed2* loc2.1* delta3* speed3* loc3.1* 
  #   0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 1.2e-05 0.0e+00 0.0e+00 
  
  
  
  
  
  ##==============================================================##
  
  
  
  
  #### ========  SPECIFICATION - ORDER 4    ======== ####
  
  ##  Method: Estimate the sub-sections of interest by changing the subset in the first line below
  ##          Log the results for each test below the code.
 
  e <- e_wbc
  # Create a TV object with just delta0, pars & st:
  ST <- 1:length(e)/length(e)
  
  # Absolute maximum: loglik = -11817, but no StdErrs
  startpars <- c(1,3,0.1, -1,3,0.4, 3,3,0.6, -3,3,0.8)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single,TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)
  parScale <- c(1,1,1,0.1,1,2,0.2,3,2,0.2,3,1,0.2)
  nDeps <- rep(1e-7,length(parScale))
  #nDeps <- c(1e-5,1e-5,1e-5,1e-9,1e-5,1e-5,1e-9,1e-5,1e-5,1e-9)  
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-9)
  
  # Close to maximum: loglik = -11823, but with some StdErrs, some NaN's
  startpars <- c(1,3,0.1, -1,3,0.4, 3,3,0.6, -3,3,0.8)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single,TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)
  #parScale <- c(1,1,1,0.1,1,2,0.2,3,2,0.2,3,1,0.2)
  parScale <- c(1,1,1,0.2,1,1,0.2,3,1,0.2,3,1,0.2)
  nDeps <- rep(1e-7,length(parScale))
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7)
  
  # Still Close to maximum: loglik = -11832, but with ALL StdErrs
  startpars <- c(1,3,0.1, -1,3,0.4, 4,3,0.6, -4,3,0.8)
  TV <- newTV(st=ST,del0=1,vecpars=startpars,shape=c(TRshape$single,TRshape$single,TRshape$single,TRshape$single),speedopt=TRspeedopt$eta)
  parScale <- c(1,1,1,0.2,1,1,0.2,3,1,0.2,3,1,0.2)
  parScale <- c(1,1,1,0.2,1,1,0.2,3,1,0.2,3,1,0.2)
  nDeps <- rep(2e-5,length(parScale))
  TV$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-6)
  
  TV <- EstimateTV(e,TV)  # ,calcHess = TRUE,verbose = TRUE)
  # Estimation Outputs:
  TV$Estimated$value
  plot(TV$Estimated$condvars,type = 'l',main=paste0("Loglik: ",round(TV$Estimated$value)))
  
  TV$Estimated$delta0
  TV$Estimated$pars
 
  ## Now calculate the standard errors for the model  ##
  TV <- calcStderr_TV(e,TV)
  TV$Estimated$stderr
  
  TV <- calcParamStats_TV(TV)
  
  TV$Estimated$parsVector
  TV$Estimated$stderr
  TV$Estimated$tStat
  TV$Estimated$PValues
  
  saveRDS(TV,"Estimated_TVOrd4_WBC.RDS")
  TV <- readRDS("Estimated_TVOrd4_WBC.RDS")
  View(TV)
  
  
  if (!exists("refData")) refData <- readRDS("RefData_withGarch_WBC_2018.RDS")
  fName <- "TestStatDist_WBC.RDS"
  # Calculate the Reference test statistics:
  refTest <- list()
  refTest$LMTR2 <- myTest.TV.noGARCH.TR2(e,TV,TESTorder$H0)
  refTest$LMRobust <- myTest.TV.noGARCH.robust(e,TV,TESTorder$H0)
  TEST <- GenTestStatDist(e,TV,refData,refTest,saveas = fName,numloops = 1600)
  
  # Confirm we managed to execute more than 1000 samples for our distribution:
  TEST$FailCount_TR2
  TEST$FailCount_Robust
  

  
  #### ========  SPECIFICATION - ORDER 4  - Results  ======== ####
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 4 ---- obs 1:7027            ##
  ##
  ## Result: Pvals=[TR2: %, Robust: %], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -11816.54 - Can't calc StdErr
  ## delta0: 2.502251
  ## Pars: 
  ##             [,1]       [,2]       [,3]        [,4]
  ##  delta -0.77819918 -0.9344535 19.9508937 -19.2324422
  ##  speed  7.00000000  5.7393151  4.7095049   4.2376019
  ##  loc1   0.09372077  0.4272500  0.6041444   0.6268639
  ##  loc2           NA         NA         NA          NA
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##==============================================================##
  
  
  ##==============================================================##
  ##  TV - AUS_4Banks\WBC ---- ORDER 4 ---- obs 1:7027            ##
  ##
  ## Result: Pvals=[TR2: %, Robust: %], using
  ## TestStat_ProbDist_WBC Logliklihood Value:  -11832.14 - Can calc StdErr
  ## delta0: 2.8424
  ## Pars: 
  ##             [,1]       [,2]       [,3]        [,4]
  ##  delta -1.3432191 -2.8119254 18.5712063 -16.3580659
  ##  speed  0.6237312  1.9862233  4.6116039   5.2115295
  ##  loc1   0.1017319  0.5735573  0.6088656   0.6326302
  ##  loc2          NA         NA         NA          NA
  ## TR2 - Cannot Reject the null: H0=Model is sufficient
  ## Robust - Cannot Reject the null: H0=Model is sufficient
  ##
  ## TV$Estimated$stderr
  ##  delta0   delta1   speed1   loc1.1   delta2   speed2   loc2.1   delta3   speed3   loc3.1   delta4   speed4   loc4.1 
  ##  1.551266 2.841646 1.790091 0.448503 0.613840 0.199490 0.025026 6.776948 0.120729 0.005850 6.796891 0.403653 0.003182 
  #
  ## TV$Estimated$tStat
  ##  delta0     delta1     speed1     loc1.1     delta2     speed2     loc2.1     delta3     speed3     loc3.1     delta4     speed4     loc4.1 
  ##  1.832310  -0.472690   0.348435   0.226825  -4.580876   9.956505  22.918884   2.740349  38.198005 104.085694  -2.406698  12.910918 198.828676 
  #
  ## TV$Estimated$PValues
  ##  delta0   delta1   speed1   loc1.1  delta2*  speed2*  loc2.1*  delta3*  speed3*  loc3.1*  delta4*  speed4*  loc4.1* 
  ##  0.063600 0.604627 0.691148 0.779538 0.000004 0.000000 0.000000 0.005845 0.000000 0.000000 0.015317 0.000000 0.000000 
  #
  ##==============================================================##
  
  
  #### ====================== Conclusion ===================== ####
  
  ## Recommend the Order 3 model
 
  
  ##==============================================================##
  ##                            THE END
  ##==============================================================##
  
  