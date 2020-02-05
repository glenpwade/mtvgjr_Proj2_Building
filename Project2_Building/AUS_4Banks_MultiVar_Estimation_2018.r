##==========================================================================================##
##  Note:  The Multivar Specification has not been rigorously tested at this stage, so...
##         This is a 'quick n dirty' estimation just to get indicative results.
##         We will complete the full specification & estimation later.  23-Feb-2018
##==========================================================================================##

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


####  ====================== Initialisation ========================  ####
rm(list=ls())
gc(T)

if (TRUE){
  library(graphics)
  library(foreach)
  library(RevoUtilsMath)

  #setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC - GOOGLE DRIVE
  #setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's laptop - GOOGLE DRIVE
  setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Glens Home PC
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v7.r")
  source(functionsFile)

}
##====================== Initialisation ========================##

##====================== Data Setup ============================##
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
    mydata <- readRDS("AUS_4Banks-ReturnsData-2018.RDS")
    dates <- readRDS("AUS_4Banks-Dates-2018.RDS")
    
    e_anz <- mydata$e_anz
    e_cba <- mydata$e_cba
    e_nab <- mydata$e_nab
    e_wbc <- mydata$e_wbc
    
    rm(mydata)
    
    ANZ <- readRDS("Estimated_TVOrd3_ANZ.RDS")  
    CBA <- readRDS("Estimated_TVOrd3_CBA.RDS")  
    NAB <- readRDS("Estimated_TVOrd3_NAB.RDS")  
    WBC <- readRDS("Estimated_TVOrd3_WBC.RDS")  
    
    # Set up the OptimControls to handle Variance Targetting:
    tmp <- ANZ$optimcontrol
    ANZ$optimcontrol_long <- tmp
    tmp$parscale <- tail(tmp$parscale, -1)
    tmp$ndeps<- tail(tmp$ndeps, -1)
    ANZ$optimcontrol_short <- tmp
    #
    tmp <- CBA$optimcontrol
    CBA$optimcontrol_long <- tmp
    tmp$parscale <- tail(tmp$parscale, -1)
    tmp$ndeps<- tail(tmp$ndeps, -1)
    CBA$optimcontrol_short <- tmp
    #
    tmp <- NAB$optimcontrol
    NAB$optimcontrol_long <- tmp
    tmp$parscale <- tail(tmp$parscale, -1)
    tmp$ndeps<- tail(tmp$ndeps, -1)
    NAB$optimcontrol_short <- tmp
    #
    tmp <- WBC$optimcontrol
    WBC$optimcontrol_long <- tmp
    tmp$parscale <- tail(tmp$parscale, -1)
    tmp$ndeps<- tail(tmp$ndeps, -1)
    WBC$optimcontrol_short <- tmp

    nTV <- list()
    nTV[["ANZ"]] <- ANZ
    nTV[["CBA"]] <- CBA
    nTV[["NAB"]] <- NAB
    nTV[["WBC"]] <- WBC
    
    saveRDS(nTV,"AUS_4Banks_nTV_2018.RDS")    
    View(nTV)
    
  }

##====================== Data Load =============================##

####================== Plot the Data ===================== ####
DoThis <- FALSE
if (DoThis) {
  ptitle <- "ANZ de-meaned Returns"
  plot(e_anz,type="l",main = ptitle)  # Index on x-axis
  plot(dates,e_anz,"l",main = ptitle) # Dates on x-axis
  
  ptitle <- "CBA de-meaned Returns"
  plot(dates,e_cba,"l",main = ptitle)
  plot(e_cba,type="l",main = ptitle)
  
  ptitle <- "NAB de-meaned Returns"
  plot(e_nab,type="l",main = ptitle)
  plot(dates,e_nab,"l",main = ptitle)
  
  ptitle <- "WBC de-meaned Returns"
  plot(e_wbc,type="l",main = ptitle)
  plot(dates,e_wbc,"l",main = ptitle)
  

}
##====================== Plot the Data =========================##


##==============================================================##
#### Step 1:  Load the 4 univariate TV objects              ####
##==============================================================##

##==============================================================##
## Univariate TV Specifications are saved in AUS_4Banks_nTV.RDS:

## See the "AUS_4Banks_UniVar-Specification_XXX_Ver2.r" files 
##         (where XXX = ANZ,CBA,NAB,WBC)
## for the full development of the univariate specifications.

nTV <- readRDS("AUS_4Banks_nTV_2018.RDS")

TV_anz <- nTV$ANZ
TV_cba <- nTV$CBA
TV_nab <- nTV$NAB
TV_wbc <- nTV$WBC

N <- length(nTV)

##==============================================================##
#### Step 2:  Create the 4 univariate GARCH objects           ####
##==============================================================##

# Note: This step is necessary for the multivariate estimation,
#       because it helps provide us with suitable starting
#       parameters and optimcontrol values


## ANZ: Garch Univar Estimateion  ####

g <- TV_anz$Estimated$condvars
w <- e_anz/sqrt(g)
# Setup GARCH object with default type & starting parameters
GARCH <- newGARCH(GARCHtype$GJR,c(0.05,0.05,0.85,0.05))

# Estimate the univar series:
#GARCH$optimcontrol <- list(fnscale = -1,parscale=c(1,1,3,1),ndeps=c(1e-4,1e-4,1e-9,1e-4),reltol=1e-8)
GARCH <- EstimateGARCH(w,GARCH,calcHess = TRUE)
h <- GARCH$Estimated$condvars

# Look at the plots:
plot(e_anz,type="l",main="ANZ Returns")
abline(v=seq(1,NROW(e_anz),by = 100),lwd=0.25,col="lightgrey")
lines(e_anz,type="l")

plot(g,type="l",main="ANZ TV(g)")
abline(v=seq(1,NROW(e_anz),by = 100),lwd=0.25,col="lightgrey")
lines(g,type="l")

plot(w,type="l",main="ANZ Returns / sqrt(g)")
abline(v=seq(1,NROW(e_anz),by = 100),lwd=0.25,col="lightgrey")
lines(w,type="l")

plot(h,type="l",main = "ANZ Garch(h)")
abline(v=seq(1,NROW(e_anz),by = 100),lwd=0.25,col="lightgrey")
lines(h,type="l")

## Now compare with GARCH Estimate for demeaned returns:
GARCH0 <- EstimateGARCH(e_anz,GARCH,calcHess = TRUE)
h0 <- GARCH0$Estimated$condvars
plot(h0,type="l",main = "ANZ Garch(h0) - with TV")
abline(v=seq(1,NROW(e_anz),by = 100),lwd=0.25,col="lightgrey")
lines(h0,type="l")

ptitle <- "GARCH_Univar_Estimate_ANZ - with TV"
cat("\n", ptitle, " - Value: ", GARCH0$Estimated$value, "\nPars:",GARCH0$Estimated$pars, "\n\n")

ptitle <- "GARCH_Univar_Estimate_ANZ - TV removed"
cat("\n", ptitle, " - Value: ", GARCH$Estimated$value, "\nPars:",GARCH$Estimated$pars, "\n\n")

## Calculate StdErrors:
GARCH0 <- calcStderr_GARCH(e_anz,GARCH0)
GARCH0$Estimated$stderr
GARCH <- calcStderr_GARCH(w,GARCH)
GARCH$Estimated$stderr

GARCH0 = calcParamStats_GARCH(GARCH0)
GARCH = calcParamStats_GARCH(GARCH)

GARCH0$Estimated$parsVector
GARCH0$Estimated$stderr
GARCH0$Estimated$PValues

GARCH$Estimated$parsVector
GARCH$Estimated$stderr
GARCH$Estimated$PValues

GARCH_anz <- GARCH


## CBA: Garch Univar Estimateion  ####
g <- TV_cba$Estimated$condvars
w <- e_cba/sqrt(g)
# Setup GARCH object with default type & starting parameters
GARCH <- newGARCH(GARCHtype$GJR,c(0.05,0.05,0.85,0.05))

# Estimate the univar series:
#GARCH$optimcontrol <- list(fnscale = -1,parscale=c(1,1,3,1),ndeps=c(1e-4,1e-4,1e-9,1e-4),reltol=1e-8)
GARCH <- EstimateGARCH(w,GARCH,calcHess = TRUE)
h <- GARCH$Estimated$condvars

# Look at the plots:
plot(e_cba,type="l",main="CBA Returns")
abline(v=seq(1,NROW(e_cba),by = 100),lwd=0.25,col="lightgrey")
lines(e_cba,type="l")

plot(g,type="l",main="CBA TV(g)")
abline(v=seq(1,NROW(e_cba),by = 100),lwd=0.25,col="lightgrey")
lines(g,type="l")

plot(w,type="l",main="CBA Returns / sqrt(g)")
abline(v=seq(1,NROW(e_cba),by = 100),lwd=0.25,col="lightgrey")
lines(w,type="l")

plot(h,type="l",main = "CBA Garch(h)")
abline(v=seq(1,NROW(e_cba),by = 100),lwd=0.25,col="lightgrey")
lines(h,type="l")

## Now compare with GARCH Estimate for demeaned returns:
GARCH0 <- EstimateGARCH(e_cba,GARCH,calcHess = TRUE)
h0 <- GARCH0$Estimated$condvars
plot(h0,type="l",main = "CBA Garch(h0) - with TV")
abline(v=seq(1,NROW(e_cba),by = 100),lwd=0.25,col="lightgrey")
lines(h0,type="l")

ptitle <- "GARCH_Univar_Estimate_CBA - with TV"
cat("\n", ptitle, " - Value: ", GARCH0$Estimated$value, "\nPars:",GARCH0$Estimated$pars, "\n\n")

ptitle <- "GARCH_Univar_Estimate_CBA - TV removed"
cat("\n", ptitle, " - Value: ", GARCH$Estimated$value, "\nPars:",GARCH$Estimated$pars, "\n\n")

## Calculate StdErrors:
GARCH0 <- calcStderr_GARCH(e_cba,GARCH0)
GARCH0$Estimated$stderr
GARCH <- calcStderr_GARCH(w,GARCH)
GARCH$Estimated$stderr

GARCH0 = calcParamStats_GARCH(GARCH0)
GARCH = calcParamStats_GARCH(GARCH)

GARCH0$Estimated$parsVector
GARCH0$Estimated$stderr
GARCH0$Estimated$PValues

GARCH$Estimated$parsVector
GARCH$Estimated$stderr
GARCH$Estimated$PValues

GARCH_cba <- GARCH

## NAB: Garch Univar Estimateion  ####
g <- TV_nab$Estimated$condvars
w <- e_nab/sqrt(g)
# Setup GARCH object with default type & starting parameters
GARCH <- newGARCH(GARCHtype$GJR,c(0.05,0.05,0.85,0.05))

# Estimate the univar series:
GARCH$optimcontrol <- list(fnscale = -1,parscale=c(1,1,3,1),ndeps=c(1e-4,1e-4,1e-7,1e-4),reltol=1e-8)
GARCH <- EstimateGARCH(w,GARCH,calcHess = TRUE)
h <- GARCH$Estimated$condvars

# Look at the plots:
plot(e_nab,type="l",main="NAB Returns")
abline(v=seq(1,NROW(e_nab),by = 100),lwd=0.25,col="lightgrey")
lines(e_nab,type="l")

plot(g,type="l",main="NAB TV(g)")
abline(v=seq(1,NROW(e_nab),by = 100),lwd=0.25,col="lightgrey")
lines(g,type="l")

plot(w,type="l",main="NAB Returns / sqrt(g)")
abline(v=seq(1,NROW(e_nab),by = 100),lwd=0.25,col="lightgrey")
lines(w,type="l")

plot(h,type="l",main = "NAB Garch(h)")
abline(v=seq(1,NROW(e_nab),by = 100),lwd=0.25,col="lightgrey")
lines(h,type="l")

## Now compare with GARCH Estimate for demeaned returns:
GARCH0 <- EstimateGARCH(e_nab,GARCH,calcHess = TRUE)
h0 <- GARCH0$Estimated$condvars
plot(h0,type="l",main = "NAB Garch(h0) - with TV")
abline(v=seq(1,NROW(e_nab),by = 100),lwd=0.25,col="lightgrey")
lines(h0,type="l")

ptitle <- "GARCH_Univar_Estimate_NAB - with TV"
cat("\n", ptitle, " - Value: ", GARCH0$Estimated$value, "\nPars:",GARCH0$Estimated$pars, "\n\n")

ptitle <- "GARCH_Univar_Estimate_NAB - TV removed"
cat("\n", ptitle, " - Value: ", GARCH$Estimated$value, "\nPars:",GARCH$Estimated$pars, "\n\n")

## Calculate StdErrors:
GARCH0 <- calcStderr_GARCH(e_nab,GARCH0)
GARCH0$Estimated$stderr
GARCH <- calcStderr_GARCH(w,GARCH)
GARCH$Estimated$stderr

GARCH0 = calcParamStats_GARCH(GARCH0)
GARCH = calcParamStats_GARCH(GARCH)

GARCH0$Estimated$parsVector
GARCH0$Estimated$stderr
GARCH0$Estimated$PValues

GARCH$Estimated$parsVector
GARCH$Estimated$stderr
GARCH$Estimated$PValues

GARCH_nab <- GARCH


## WBC:Garch Univar Estimateion  ####
g <- TV_wbc$Estimated$condvars
w <- e_wbc/sqrt(g)
# Setup GARCH object with default type & starting parameters
GARCH <- newGARCH(GARCHtype$GJR,c(0.05,0.05,0.85,0.05))

# Estimate the univar series:
GARCH$optimcontrol <- list(fnscale = -1,parscale=c(1,1,3,1),ndeps=c(1e-4,1e-4,1e-7,1e-4),reltol=1e-8)
GARCH <- EstimateGARCH(w,GARCH,calcHess = TRUE)
h <- GARCH$Estimated$condvars

# Look at the plots:
plot(e_wbc,type="l",main="WBC Returns")
abline(v=seq(1,NROW(e_wbc),by = 100),lwd=0.25,col="lightgrey")
lines(e_wbc,type="l")

plot(g,type="l",main="WBC TV(g)")
abline(v=seq(1,NROW(e_wbc),by = 100),lwd=0.25,col="lightgrey")
lines(g,type="l")

plot(w,type="l",main="WBC Returns / sqrt(g)")
abline(v=seq(1,NROW(e_wbc),by = 100),lwd=0.25,col="lightgrey")
lines(w,type="l")

plot(h,type="l",main = "WBC Garch(h)")
abline(v=seq(1,NROW(e_wbc),by = 100),lwd=0.25,col="lightgrey")
lines(h,type="l")

## Now compare with GARCH Estimate for demeaned returns:
GARCH0 <- EstimateGARCH(e_wbc,GARCH,calcHess = TRUE)
h0 <- GARCH0$Estimated$condvars
plot(h0,type="l",main = "WBC Garch(h0) - with TV")
abline(v=seq(1,NROW(e_wbc),by = 100),lwd=0.25,col="lightgrey")
lines(h0,type="l")

ptitle <- "GARCH_Univar_Estimate_WBC - with TV"
cat("\n", ptitle, " - Value: ", GARCH0$Estimated$value, "\nPars:",GARCH0$Estimated$pars, "\n\n")

ptitle <- "GARCH_Univar_Estimate_WBC - TV removed"
cat("\n", ptitle, " - Value: ", GARCH$Estimated$value, "\nPars:",GARCH$Estimated$pars, "\n\n")

## Calculate StdErrors:
GARCH0 <- calcStderr_GARCH(e_wbc,GARCH0)
GARCH0$Estimated$stderr
GARCH <- calcStderr_GARCH(w,GARCH)
GARCH$Estimated$stderr

GARCH0 = calcParamStats_GARCH(GARCH0)
GARCH = calcParamStats_GARCH(GARCH)

GARCH0$Estimated$parsVector
GARCH0$Estimated$stderr
GARCH0$Estimated$PValues

GARCH$Estimated$parsVector
GARCH$Estimated$stderr
GARCH$Estimated$PValues

GARCH_wbc <- GARCH


##==============================================================##
#### Step 3:  Create the multivariate objects for TV & GARCH  ####
##==============================================================##

## nTV was loaded above & used to create nGARCH

nGARCH <- list()
nGARCH$ANZ <- GARCH_anz
nGARCH$CBA <- GARCH_cba
nGARCH$NAB <- GARCH_nab
nGARCH$WBC <- GARCH_wbc
saveRDS(nGARCH,"AUS_4Banks_nGARCH_2018.RDS")

View(nGARCH)

##==============================================================##
#### Step 4:  Create the STCC object                          ####
##==============================================================##

e <- cbind(e_anz,e_cba,e_nab,e_wbc)
Tobs <- NROW(e)    

STCC <- newSTCC(Tobs)
# Set the starting parameters
STCC$P1 <- matrix(0.2,N,N)  # Const Corr = 0.2
diag(STCC$P1) <- 1
STCC$P2 <- matrix(0.8,N,N)  # Const Corr = 0.8
diag(STCC$P2) <- 1
# Estimate the model
nDeps <- c(rep(1e-4,12),1e-5,1e-4)
#parScale <- c(rep(2,12), c(3,2))
parScale <- rep(1,14)
STCC$optimcontrol <- list(fnscale = -1,parscale=parScale,ndeps=nDeps,reltol=1e-7)
STCC <- EstimateSTCC(e,STCC,calcHess = TRUE,verbose = TRUE)

# prev best value: -16628.88

saveRDS(STCC,"AUS_4Banks_STCC_2018.RDS")
View(STCC)

STCC <- readRDS("AUS_4Banks_STCC_2018.RDS")

##===================================================================##
#### Step 5:  Get initial Multivar Estimate of STCC, with h(t)=1   ####
##===================================================================##

# Method:
#
# 1: Get initial estimate of correlation, assuming No Garch, i.e. h(t)=1, for all N
# Initialise STCC Estimate (h(tn) = 1, all n = 1..N) for the multivariate loop:
#
# We use nGarch..h(t) = 1, nTV..g(t) = Estimated$condvars, and STCC = Starting values
#

# Calculate how many params are in the full model
TVPars.Count <- 1  # delta0
for (n in 1:NROW(nTV)) TVPars.Count <- TVPars.Count + length(nTV[[n]]$Estimated$parsVector)
GarchPars.Count <- 0
for (n in 1:length(nGARCH)) GarchPars.Count <- GarchPars.Count + length(nGARCH[[n]]$Estimated$parsVector)
STCCPars.Count <- length(STCC$Estimated$parsVector)
#
ModelPars.Count <- TVPars.Count + GarchPars.Count + STCCPars.Count

# Set the variance targetting mode:
V_Target <- FALSE

# Get the univariate objects:
nTV <- readRDS("AUS_4Banks_nTV_2018.RDS")
nGARCH <- readRDS("AUS_4Banks_nGARCH_2018.RDS")
nGARCH_ht1 <- nGARCH
for (n in 1:length(nGARCH_ht1)) nGARCH_ht1[[n]]$Estimated$condvars <- rep(1,Tobs)

STCC <- newSTCC(Tobs)
# Set the starting parameters
STCC$P1 <- matrix(0.2,N,N)  # Const Corr = 0.2
diag(STCC$P1) <- 1
STCC$P2 <- matrix(0.8,N,N)  # Const Corr = 0.8
diag(STCC$P2) <- 1

# Setup the multivariate Return object:
mvRTN <- list()     
mvRTN$ngarch <- nGARCH_ht1
mvRTN$ntv <- nTV
mvRTN$stcc <- STCC
mvRTN$value <- -1e10
mvRTN$m_pars <- rep(0.5,ModelPars.Count)  # Initialise the model params vector
# Setup an identical object to store the previous state:
mvRTN_prev <- mvRTN

#Build the matrix of data:
e <- cbind(e_anz,e_cba,e_nab,e_wbc)
Tobs <- NROW(e)
N <- NCOL(e)

# Now do the Estimation:
timestamp()
mvRTN <- Estimate_MTVGARCH(e,nTV,nGARCH_ht1,STCC,focus="STCC",var_target=V_Target)
timestamp()

saveRDS(mvRTN$stcc,"AUS_4Banks_STCC_ht1_2018.RDS")
View(mvRTN)

identical(STCC$Estimated$condcorrs ,mvRTN$stcc$Estimated$condcorrs)

mean(STCC$Estimated$condcorrs-mvRTN$stcc$Estimated$condcorrs)

##==============================================================##
####  Step 6: Execute multivar estimation loop  ####
##==============================================================##

# Method:
#
# 1: Get initial estimate of correlation, assuming No Garch, i.e. h(t)=1, for all N
# 2: Now estimate the nGarch, using the STCC estimates from above
# 3: Now estimate the nTV, using the nGARCH estimates from above
# 4: Estimate STCC using the nGARCH object from the previous estimate 
##  Repeat steps 2 - 4 until optimisation is complete!

# Load up the univariate objects:
nTV <- readRDS("AUS_4Banks_nTV_2018.RDS")
nGARCH <- readRDS("AUS_4Banks_nGARCH_2018.RDS")
STCC <- readRDS("AUS_4Banks_STCC_ht1_2018.RDS")

nGARCH_ht1 <- nGARCH
for (n in 1:length(nGARCH_ht1)) nGARCH_ht1[[n]]$Estimated$condvars <- rep(1,Tobs)


# Setup the multivariate Return object:
mvRTN <- list()     
mvRTN$ngarch <- nGARCH_ht1
mvRTN$ntv <- nTV
mvRTN$stcc <- STCC
mvRTN$value <- -1e10
mvRTN$m_pars <- rep(0.5,ModelPars.Count)  # Initialise the model params vector
# Setup an identical object to store the previous state:
mvRTN_prev <- mvRTN

#Build the matrix of data:
e <- cbind(e_anz,e_cba,e_nab,e_wbc)
Tobs <- NROW(e)
N <- NCOL(e)


# Set the maximum number of loops
loops <- 21
# Setup a list to track LL & params for each loop:
m_results <- list()
listIdx <- 0
for (n in 1:(3*loops)) m_results[[n]] <- list()

# Setup tolerances to exit loop if required.
value_tol <- 1e-8
param_tol <- 1e-6

# Set the variance targetting mode:
V_Target <- TRUE

timestamp()
for (i in 1:loops) {
  
  # 1: Estimate the Garch, using the STCC estimates from above
  nTV <- mvRTN$ntv
  nGARCH <- mvRTN$ngarch
  STCC <- mvRTN$stcc
  mvRTN <- EstimateM_TVGARCH(e,nTV,nGARCH,STCC,focus="GARCH",var_target=V_Target)

  listIdx <- (i-1)*3 + 1
  m_results[[listIdx]]$value <- mvRTN$value
  m_results[[listIdx]]$m_pars <- mvRTN$m_pars
  g_n <- matrix(0,Tobs,N)
  for(n in 1:N) g_n[,n] <- mvRTN$ntv[[n]]$Estimated$condvars
  m_results[[listIdx]]$g_n <- g_n
  h_n <-  matrix(0,Tobs,N)
  for(n in 1:N) h_n[,n] <- mvRTN$ngarch[[n]]$Estimated$condvars
  m_results[[listIdx]]$h_n <- h_n
  # Cond. Corr's  
  m_results[[listIdx]]$c_n <- mvRTN$stcc$Estimated$condcorrs
  cat("\n",date(),"\nGarch done: LL Value:",mvRTN$value)
  
  # 2: Now estimate the TV, using the GARCH estimates from above
  nTV <- mvRTN$ntv
  nGARCH <- mvRTN$ngarch
  STCC <- mvRTN$stcc
  mvRTN <- EstimateM_TVGARCH(e,nTV,nGARCH,STCC,focus="TV",var_target=V_Target)

  listIdx <- (i-1)*3 + 2
  m_results[[listIdx]]$value <- mvRTN$value
  m_results[[listIdx]]$m_pars <- mvRTN$m_pars
  g_n <- matrix(0,Tobs,N)
  for(n in 1:N) g_n[,n] <- mvRTN$ntv[[n]]$Estimated$condvars
  m_results[[listIdx]]$g_n <- g_n
  h_n <-  matrix(0,Tobs,N)
  for(n in 1:N) h_n[,n] <- mvRTN$ngarch[[n]]$Estimated$condvars
  m_results[[listIdx]]$h_n <- h_n
  # Cond. Corr's  
  m_results[[listIdx]]$c_n <- mvRTN$stcc$Estimated$condcorrs
  cat("\n",date(),"\nTV done: LL Value:",mvRTN$value)
  
  # 3: Now estimate the Correlation, using the TV & GARCH estimates from above
  nTV <- mvRTN$ntv
  nGARCH <- mvRTN$ngarch
  STCC <- mvRTN$stcc
  mvRTN <- EstimateM_TVGARCH(e,nTV,nGARCH,STCC,focus="STCC",var_target=V_Target)

  listIdx <- (i-1)*3 + 3
  m_results[[listIdx]]$value <- mvRTN$value
  m_results[[listIdx]]$m_pars <- mvRTN$m_pars
  g_n <- matrix(0,Tobs,N)
  for(n in 1:N) g_n[,n] <- mvRTN$ntv[[n]]$Estimated$condvars
  m_results[[listIdx]]$g_n <- g_n
  h_n <-  matrix(0,Tobs,N)
  for(n in 1:N) h_n[,n] <- mvRTN$ngarch[[n]]$Estimated$condvars
  m_results[[listIdx]]$h_n <- h_n
  # Cond. Corr's  
  m_results[[listIdx]]$c_n <- mvRTN$stcc$Estimated$condcorrs
  cat("\n",date(),"\nSTCC done: LL Value:",mvRTN$value)

  # 4. Update & save results from this loop!
  ###  Check if things have changed!!!  ###
  cat("\n\nLL Increasing? ",isTRUE(mvRTN$value > mvRTN_prev$value), all.equal.numeric(mvRTN$value,mvRTN_prev$value,tolerance = value_tol))
  cat("\nChange in params? ", all.equal.numeric(mvRTN$m_pars,mvRTN_prev$m_pars,tolerance = param_tol))
  cat ("\nEnd of Loop number:", i,"\n")
  # Save to file!!
  saveRDS(m_results,"TVGJR_MultivarEstimation_2018_1.RDS")

  
  # Exit the loop if the parameters have stopped moving
  try (  if (isTRUE(all.equal.numeric(mvRTN$m_pars,mvRTN_prev$m_pars,tolerance = param_tol))) {
    print("Seems we could have stopped here!")
  } )
  
  mvRTN_prev <- mvRTN
  #
} # End of for (i in 1:loops)
timestamp()

saveRDS(mvRTN,"TVGJR_MultivarEstimation_2018_3.RDS")


##==============================================================##
##                            THE END
##==============================================================##


## ---------------  Analysis  ------------- ##

(m_results[[6]]$m_pars)

all.equal.numeric(m_results[[4]]$m_pars,m_results[[5]]$m_pars)
identical(m_results[[5]]$m_pars,m_results[[7]]$m_pars)

round(m_results[[5]]$m_pars,5)

# -- #
var_target <- TRUE
tmp <- m_results[[5]]$m_pars
tmpTV <- tmp[1:49]
tmpG <- tmp[50:65]
tmpST <- tail(tmp,14)

ntv <- mvRTN$ntv
ngarch <- mvRTN$ngarch
stcc <- mvRTN$stcc

focus <- "TV"
rtnTV <- myLogLik.multivar.TVGARCHSTCC(tmpTV,e,ntv,ngarch,stcc,focus,var_target,return_ll = FALSE)
focus <- "GARCH"
rtnGARCH <- myLogLik.multivar.TVGARCHSTCC(tmpG,e,ntv,ngarch,stcc,focus,var_target,return_ll = FALSE)
focus <- "STCC"
rtnSTCC <- myLogLik.multivar.TVGARCHSTCC(tmpST,e,ntv,ngarch,stcc,focus,var_target,return_ll = FALSE)



## ---------------  End: Analysis  ------------- ##


plot(STCC$Estimated$condcorrs[,1],type = 'l',ylim=c(0.45,0.85),panel.first = grid())
lines(STCC$Estimated$condcorrs[,2],type = 'l',col="red")
lines(STCC$Estimated$condcorrs[,3],type = 'l',col="green")
lines(STCC$Estimated$condcorrs[,4],type = 'l',col="blue")
lines(STCC$Estimated$condcorrs[,5],type = 'l',col="darkgrey")
lines(STCC$Estimated$condcorrs[,6],type = 'l',col="yellow",lwd=2)

saveRDS(mvRTN,"MTVGJR_Estimated_vartarget_OFF_2018.RDS")
mvRTN <- readRDS("MTVGJR_Estimated_vartarget_OFF_2018.RDS")

plot(e_anz,type='l')

## Run a rolling window over the data - calc local correlation ##

calculate_local_corr <- function(e,localvarwindow=500) {

  Tobs <- NROW(e)
  localCorr <- matrix(NA,Tobs,6)
  
  for (t in 1:Tobs) {
    if (t <= (localvarwindow)) {
      localCorr[t,] <- as.vector(myVecl(cor(e[1:localvarwindow,])))
    } else if (t > localvarwindow && t < (Tobs-(localvarwindow+1))) {
      localCorr[t,] <- as.vector(myVecl(cor(e[t:(t+localvarwindow),])))
    } else localCorr[t,] <- as.vector(myVecl(cor(e[(Tobs-localvarwindow):Tobs,])))
  }
  
  # Return:
  localCorr
}

e <- cbind(e_anz,e_cba,e_nab,e_wbc)

lcor <- calculate_local_corr(e,250)

plot(mvRTN$stcc$Estimated$condcorrs[,1],col="blue",ylab="",xlab="",lwd=2,type = 'l',ylim=c(0,1),panel.first = grid())
#plot(lcor[,1],type = 'l')
lines(lcor[,1],type = 'l')
lines(lcor[,2],type = 'l')
lines(lcor[,3],type = 'l')
lines(lcor[,4],type = 'l')
lines(lcor[,5],type = 'l')
lines(lcor[,6],type = 'l')
lines((nTV$ANZ$Estimated$condvars/12),type = 'l',col="red")
legend("topleft",inset = 0.05,legend = "Window: 250", box.col = "white")

## -----------------------------------------------   ##


####  ======  Calculate Standard Errors  ======  ####





