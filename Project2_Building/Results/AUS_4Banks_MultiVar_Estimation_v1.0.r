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
  library(doParallel)

  #setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC - GOOGLE DRIVE
  setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's laptop - GOOGLE DRIVE
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v6.r")
  source(functionsFile)

}
##====================== Initialisation ========================##

##====================== Data Setup ============================##
DoThis <- FALSE
if (DoThis){
  mydata <- read.csv("Bank-Returns.csv",header=TRUE)
  dates <- as.Date(mydata$Date, format = "%d/%m/%Y")

  #Create scaled-up individtal data series vectors (Percentage Returns):
  mydata$e_anz <- mydata$anz * 100
  mydata$e_cba <- mydata$cba * 100
  mydata$e_nab <- mydata$nab * 100
  mydata$e_wbc <- mydata$wbc * 100
  
  #De-mean the returns data:
  mydata$e_anz <- mydata$e_anz - mean(mydata$e_anz)
  mydata$e_cba <- mydata$e_wbc - mean(mydata$e_cba)
  mydata$e_nab <- mydata$e_nab - mean(mydata$e_nab)
  mydata$e_wbc <- mydata$e_wbc - mean(mydata$e_wbc)
  
  # save data  
  saveRDS(mydata,"AUS_4Banks-ReturnsData.RDS")
  saveRDS(dates,"AUS_4Banks-Dates.RDS")
}
##====================== Data Setup ============================##

##====================== Data Load =============================##
if (TRUE){
  # Read AUS_4Banks data from saved file
  dates <- readRDS("AUS_4Banks-Dates.RDS")
  mydata <- readRDS("AUS_4Banks-ReturnsData.RDS")
  e_anz <- mydata$e_anz
  e_cba <- mydata$e_cba
  e_nab <- mydata$e_nab
  e_wbc <- mydata$e_wbc
}
##====================== Data Load =============================##

##====================== Plot the Data =========================##
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

nTV <- readRDS("AUS_4Banks_nTV.RDS")

TV_anz <- nTV$ANZ
TV_cba <- nTV$CBA
TV_nab <- nTV$NAB
TV_wbc <- nTV$WBC

N <- length(nTV)

##==============================================================##
#### Step 2:  Create the 4 univariate GARCH objects           ####
##==============================================================##

# Note: This step is not necessary for the multivariate estimation,
#       but it may be useful in providing us with suitable starting
#       parameters...


# ANZ:

e <- e_anz/sqrt(TV_anz$Estimated$condvars)
# Setup GARCH object with starting parameters
GARCH <- list()
GARCH$type <- GARCHtype$GJR
GARCH$pars <- c(0.05,0.05,0.85,0.05)  #Starting params
GARCH$var_target <- FALSE   # Omega (param #1) is free, not calculated based on alpha & beta
parScale <- c(3,3,1,3)
nDeps <- rep(1e-4,length(parScale))
GARCH$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
# Estimate the univar series:
GARCH <- EstimateGARCH(e,GARCH)
ptitle <- "GARCH_Univar_Estimate_ANZ"
cat("\n", ptitle, "Logliklihood Value: ", GARCH$Estimated$value, "\nPars:",GARCH$Estimated$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(GARCH$Estimated$condvars,type="l")

GARCH_anz <- GARCH


# CBA:

e <- e_cba/sqrt(TV_cba$Estimated$condvars)
# Setup GARCH object with starting parameters
GARCH <- list()
GARCH$type <- GARCHtype$GJR
GARCH$pars <- c(0.05,0.05,0.85,0.05)  #Starting params
GARCH$var_target <- FALSE   # Omega (param #1) is free, not calculated based on alpha & beta
parScale <- c(3,3,1,3)
nDeps <- rep(1e-4,length(parScale))
GARCH$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
# Estimate the univar series:
GARCH <- EstimateGARCH(e,GARCH)
ptitle <- "GARCH_Univar_Estimate_CBA"
cat("\n", ptitle, "Logliklihood Value: ", GARCH$Estimated$value, "\nPars:",GARCH$Estimated$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(GARCH$Estimated$condvars,type="l")

GARCH_cba <- GARCH


# NAB:

e <- e_nab/sqrt(TV_nab$Estimated$condvars)
# Setup GARCH object with starting parameters
GARCH <- list()
GARCH$type <- GARCHtype$GJR
GARCH$pars <- c(0.05,0.05,0.85,0.05)  #Starting params
GARCH$var_target <- FALSE   # Omega (param #1) is free, not calculated based on alpha & beta
parScale <- c(3,3,1,3)
nDeps <- rep(1e-4,length(parScale))
GARCH$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
# Estimate the univar series:
GARCH <- EstimateGARCH(e,GARCH)
ptitle <- "GARCH_Univar_Estimate_NAB"
cat("\n", ptitle, "Logliklihood Value: ", GARCH$Estimated$value, "\nPars:",GARCH$Estimated$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(GARCH$Estimated$condvars,type="l")

GARCH_nab <- GARCH


# WBC:

e <- e_wbc/sqrt(TV_wbc$Estimated$condvars)
# Setup GARCH object with starting parameters
GARCH <- list()
GARCH$type <- GARCHtype$GJR
GARCH$pars <- c(0.05,0.05,0.85,0.05)  #Starting params
GARCH$var_target <- FALSE   # Omega (param #1) is free, not calculated based on alpha & beta
parScale <- c(3,3,1,3)
nDeps <- rep(1e-4,length(parScale))
GARCH$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps = nDeps)
# Estimate the univar series:
GARCH <- EstimateGARCH(e,GARCH)
ptitle <- "GARCH_Univar_Estimate_WBC"
cat("\n", ptitle, "Logliklihood Value: ", GARCH$Estimated$value, "\nPars:",GARCH$Estimated$pars, "\n\n")
# Look at the plots:
plot(e,type="l")
plot(GARCH$Estimated$condvars,type="l")

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
saveRDS(nGARCH,"AUS_4Banks_nGARCH.RDS")


##==============================================================##
#### Step 4:  Create the STCC object                          ####
##==============================================================##

e <- cbind(e_anz,e_cba,e_nab,e_wbc)

STCC <- list()
STCC$Tobs <- NROW(e)    
STCC$type <- STCCtype$STCC
STCC$P1 <- matrix(0.2,N,N)  # Const Corr = 0.2
diag(STCC$P1) <- 1
STCC$P2 <- matrix(0.8,N,N)  # Const Corr = 0.8
diag(STCC$P2) <- 1
STCC$TRpars <- c(2,0.5)  
STCC$speedoption <- TVspeedopt$eta    
STCC$st <- 1:STCC$Tobs/STCC$Tobs 
STCC$shape <- TVshape$linear    
numCovPars <- NROW(myVecl(STCC$P1))
parScale <- rep(1,(2*numCovPars+length(STCC$TRpars)))
nDeps <- rep(1e-4,length(parScale))
STCC$optimcontrol <- list(fnscale = -1, parscale=parScale, ndeps=nDeps)

saveRDS(STCC,"AUS_4Banks_STCC.RDS")
#

##==============================================================##
#### Step 5:  Begin the multivariate estimation               ####
##==============================================================##


#Build the matrix of data:
e <- cbind(e_anz,e_cba,e_nab,e_wbc)

# Set the variance targetting mode:
V_Tartget <- TRUE

## We need to manually track the LL value now, so:
mvRTN <- NA
mvRTN_prev <- NA

# 1: Get initial estimate of correlation, assuming No Garch:

# 2: Now estimate the Garch, using the STCC estimates from above

# 3: Now estimate the TV, using the GARCH estimates from above

# 4: Estimate STCC using the nGARCH object from the previous estimate 

##  Repeat steps 2 - 4 until optimisation is complete!

M_TV_GJR <- mvRTN
saveRDS(M_TV_GJR, "M_TV_GJR.RDS")

##==============================================================##
##                            THE END
##==============================================================##


## Look at the correlations:
for (n in 1:6) {
  plot(STCC$Estimated$condcorrs[,n],lwd=1,main=paste0("Corr #",n))
}



##==============================================================##
####  Step 6: Execute multivar estimation loop  ####
##==============================================================##

# Load up the univariate objects:
nTV <- readRDS("AUS_4Banks_nTV.RDS")
nGARCH <- readRDS("AUS_4Banks_nGARCH.RDS")
STCC <- readRDS("STCC_ht1.RDS")

#Build the matrix of data:
e <- cbind(e_anz,e_cba,e_nab,e_wbc)
Tobs <- NROW(e)
N <- NCOL(e)

# Set the variance targetting mode:
V_Target <- FALSE

# Set up the multivariate Return object:
mvRTN <- list()     
mvRTN$ngarch <- nGARCH
mvRTN$ntv <- nTV
mvRTN$stcc <- STCC
mvRTN_prev <- mvRTN

# Setup a list to track LL for each loop:
loops <- 3
m_results <- list()
listIdx <- 0
for (n in 1:(3*loops)) m_results[[n]] <- list()

for (n in 1:loops) {
  
  # 1: Estimate the Garch, using the STCC estimates from above
  
  nTV <- mvRTN$ntv
  STCC <- mvRTN$stcc
  tmr <- proc.time()
  mvRTN <- EstimateM_TVGARCH(e,nTV,nGARCH,STCC,focus="GARCH",var_target=V_Target)
  proc.time() - tmr
  # 900 secs
  # Compare prev & new
  
  listIdx <- (n-1)*3 + 1
  m_results[[listIdx]]$ll_value <- mvRTN$value
  m_results[[listIdx]]$m_pars <- mvRTN$m_pars
  
  # 2: Now estimate the TV, using the GARCH estimates from above
  nGARCH <- mvRTN$ngarch
  STCC <- mvRTN$stcc
  tmr <- proc.time()
  mvRTN <- EstimateM_TVGARCH(e,nTV,nGARCH,STCC,focus="TV",var_target=V_Target)
  proc.time() - tmr
  # 140 secs
  # Compare prev & new
  
  listIdx <- (n-1)*3 + 2
  m_results[[listIdx]]$ll_value <- mvRTN$value
  m_results[[listIdx]]$m_pars <- mvRTN$m_pars
  
  # 3: Now estimate the Correlation, using the TV & GARCH estimates from above
  nTV <- mvRTN$ntv
  nGARCH <- mvRTN$ngarch
  tmr <- proc.time()
  mvRTN <- EstimateM_TVGARCH(e,nTV,nGARCH,STCC,focus="STCC",var_target=V_Target)
  proc.time() - tmr
  # Compare prev & new
  
  listIdx <- (n-1)*3 + 3
  m_results[[listIdx]]$ll_value <- mvRTN$value
  m_results[[listIdx]]$m_pars <- mvRTN$m_pars

  # 4. Update & save results from this loop!
  cat ("\nEnd of Loop number:", n)
  # Save to file!!
  saveRDS(m_results,"TVGJR_MultivarEstimation_vtF.RDS")

  ###  Check if things have changed!!!  ###
  
  identical(nTV,mvRTN$ntv)
  identical(nGARCH,mvRTN$ngarch)
  identical(STCC,mvRTN$stcc)
  identical(mvRTN$value,mvRTN_prev$value)
  identical(mvRTN$m_pars,mvRTN_prev$m_pars)
  identical(mvRTN,mvRTN_prev)
  #
}


plot(STCC$Estimated$condcorrs[,1],type = 'l',ylim=c(0.45,0.85),panel.first = grid())
lines(STCC$Estimated$condcorrs[,2],type = 'l',col="red")
lines(STCC$Estimated$condcorrs[,3],type = 'l',col="green")
lines(STCC$Estimated$condcorrs[,4],type = 'l',col="blue")
lines(STCC$Estimated$condcorrs[,5],type = 'l',col="darkgrey")
lines(STCC$Estimated$condcorrs[,6],type = 'l',col="yellow",lwd=2)

saveRDS(mvRTN,"MTVGJR_Estimated_vartarget_OFF.RDS")
mvRTN <- readRDS("MTVGJR_Estimated.RDS")


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

e <- cbind(e_anz,e_cba,e_nab,e_wbc)

ntv <- mvRTN$ntv
ngarch <- mvRTN$ngarch
stcc <- mvRTN$stcc
var_target <- TRUE
calcHess <- TRUE

## STCC Std Errors:
optimpars <- tail(mvRTN$m_pars,14)
focus <- "STCC"
hess <- optimHess(optimpars,myLogLik.multivar.TVGARCHSTCC,gr=NULL,e,ntv,ngarch,stcc,focus,var_target,return_ll=TRUE)

optimpars
se <- sqrt(-diag(solve(hess))) 


## GARCH Std Errors:
optimpars <- (mvRTN$m_pars[50:65])
focus <- "GARCH"
hess <- optimHess(optimpars,myLogLik.multivar.TVGARCHSTCC,gr=NULL,e,ntv,ngarch,stcc,focus,var_target,return_ll=TRUE)

optimpars
gse <- sqrt(-diag(solve(hess))) 

## TV Std Errors:
optimpars <- (mvRTN$m_pars[1:49])
focus <- "TV"
hess <- optimHess(optimpars,myLogLik.multivar.TVGARCHSTCC,gr=NULL,e,ntv,ngarch,stcc,focus,var_target,return_ll=TRUE)
tvse <- sqrt(-diag(solve(hess))) 
optimpars
tvse



