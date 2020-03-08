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
  library(ts)

  setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building_2020")        # Anna's work PC - GOOGLE DRIVE
  #setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building_2020")        # Anna's laptop - GOOGLE DRIVE
  
  source("clsCORR.r")
  source("clsTV.r")
  source("clsGARCH.r")

}
##====================== Initialisation ========================##

##====================== Data Setup ============================##
DoThis <- FALSE
if (DoThis){
  mydata <- read.csv("Bank-Returns-18.csv",header=TRUE)
  dates <- as.Date(mydata$date, format = "%d/%m/%Y")

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
  saveRDS(mydata,"AUS_4Banks-ReturnsData_v2.1.RDS")
  saveRDS(dates,"AUS_4Banks-Dates_v2.1.RDS")
}
##====================== Data Setup ============================##

##====================== Data Load =============================##
if (TRUE){
  # Read AUS_4Banks data from saved file
  dates <- readRDS("AUS_4Banks-Dates_v2.1.RDS")
  mydata <- readRDS("AUS_4Banks-ReturnsData_v2.1.RDS")
  e_anz <- mydata$e_anz
  e_cba <- mydata$e_cba
  e_nab <- mydata$e_nab
  e_wbc <- mydata$e_wbc
}
##====================== Data Load =============================##

##====================== Plot the Data =========================##
DoThis <- FALSE
if (DoThis) {
  ymin <- -15
  ymax <- 15
  par(mar = c(2.5, 2.5, 2, 1))
  
  ptitle <- "ANZ"
  #plot(e_anz,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))  # Index on x-axis
  plot(dates,e_anz,"l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax)) # Dates on x-axis
  # dev.copy(pdf,'ANZ_v2.pdf',width=5, height=5)
  # dev.off()
  ptitle <- "CBA"
  plot(dates,e_cba,"l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  #plot(e_cba,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  # dev.copy(pdf,'CBA_v2.pdf',width=5, height=5)
  # dev.off()
  ptitle <- "NAB"
  #plot(e_nab,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  plot(dates,e_nab,"l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  # dev.copy(pdf,'NAB_v2.pdf',width=5, height=5)
  # dev.off()
  ptitle <- "WBC"
  #plot(e_wbc,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  plot(dates,e_wbc,"l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  # dev.copy(pdf,'WBC_v2.pdf',width=5, height=5)
  # dev.off()
}
  ################# ABS(ret) and g PLOTS ########################
DoThis <- FALSE
if(DoThis){
  tmp<- readRDS("Results/TVGJR_MultivarEstimation_2018_2.RDS")
  ymin <- 0
  ymax <- 10
  par(mar = c(2.5, 2.5, 2, 1))
  
  ptitle <- "ANZ"
  g <- tmp$ntv$ANZ$Estimated$condvars
  plot(dates,abs(e_anz),type='l',col='gray80',main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  par(new = TRUE)
  plot(dates,sqrt(g),type='l',main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  # dev.copy(pdf,'ANZ_v2_g.pdf',width=5, height=5)
  # dev.off()
  ptitle <- "CBA"
  g <- tmp$ntv$CBA$Estimated$condvars
  plot(dates,abs(e_cba),type='l',col='gray80',main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  par(new = TRUE)
  plot(dates,sqrt(g),type='l',main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  # dev.copy(pdf,'CBA_v2_g.pdf',width=5, height=5)
  # dev.off()
  ptitle <- "NAB"
  g <- tmp$ntv$NAB$Estimated$condvars
  plot(dates,abs(e_nab),type='l',col='gray80',main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  par(new = TRUE)
  plot(dates,sqrt(g),type='l',main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  # dev.copy(pdf,'NAB_v2_g.pdf',width=5, height=5)
  # dev.off()
  ptitle <- "WBC"
  g <- tmp$ntv$WBC$Estimated$condvars
  plot(dates,abs(e_wbc),type='l',col='gray80',main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  par(new = TRUE)
  plot(dates,sqrt(g),type='l',main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  # dev.copy(pdf,'WBC_v2_g.pdf',width=5, height=5)
  # dev.off()
}
  ################# ACF PLOTS ########################
DoThis <- FALSE
if(DoThis){
  ymin<-0
  ymax<-0.3
  par(font.main=1)
  par(mar = c(2.5, 2.5, 2, 1))
  
  ptitle <- "ANZ"
  ser <- e_anz^2
  g <- tmp$ntv$ANZ$Estimated$condvars
  Acf(ser,lag.max=100,main="",font.main=1,xlab="", ylim=c(ymin, ymax))
  title(ptitle,line=0.5)
  #dev.copy(pdf,'ANZ_v2_acf_before.pdf',width=7, height=3)
  #dev.off()
  Acf(ser/g,lag.max=100,main="",font.main=1,xlab="", ylim=c(ymin, ymax))
  title(ptitle,line=0.5)
  #dev.copy(pdf,'ANZ_v2_acf_after.pdf',width=7, height=3)
  #dev.off()
  
  ptitle <- "CBA"
  ser <- e_cba^2
  g <- tmp$ntv$CBA$Estimated$condvars
  Acf(ser,lag.max=100,main="",font.main=1,xlab="", ylim=c(ymin, ymax))
  title(ptitle,line=0.5)
  # dev.copy(pdf,'CBA_v2_acf_before.pdf',width=7, height=3)
  # dev.off()
  Acf(ser/g,lag.max=100,main="",font.main=1,xlab="", ylim=c(ymin, ymax))
  title(ptitle,line=0.5)
  # dev.copy(pdf,'CBA_v2_acf_after.pdf',width=7, height=3)
  # dev.off()
  
  ptitle <- "NAB"
  ser <- e_nab^2
  g <- tmp$ntv$NAB$Estimated$condvars
  Acf(ser,lag.max=100,main="",font.main=1,xlab="", ylim=c(ymin, ymax))
  title(ptitle,line=0.5)
  # dev.copy(pdf,'NAB_v2_acf_before.pdf',width=7, height=3)
  # dev.off()
  Acf(ser/g,lag.max=100,main="",font.main=1,xlab="", ylim=c(ymin, ymax))
  title(ptitle,line=0.5)
  # dev.copy(pdf,'NAB_v2_acf_after.pdf',width=7, height=3)
  # dev.off()
  
  ptitle <- "WBC"
  ser <- e_wbc^2
  g <- tmp$ntv$WBC$Estimated$condvars
  Acf(ser,lag.max=100,main="",font.main=1,xlab="", ylim=c(ymin, ymax))
  title(ptitle,line=0.5)
  # dev.copy(pdf,'WBC_v2_acf_before.pdf',width=7, height=3)
  # dev.off()
  Acf(ser/g,lag.max=100,main="",font.main=1,xlab="", ylim=c(ymin, ymax))
  title(ptitle,line=0.5)
  # dev.copy(pdf,'WBC_v2_acf_after.pdf',width=7, height=3)
  # dev.off()
  
}
################# CORRELATION PLOTS ########################
DoThis <- FALSE
if(DoThis){
  tmp<- readRDS("Results/TVGJR_MultivarEstimation_2018_2.RDS")
  ymin<-0.4
  ymax<-0.9
  ptitle <- ""
  par(font.main=1)
  par(mar = c(2.5, 2.5, 2, 1))
  #View(tmp)
  cor_ANZ_CBA <- tmp$stcc$Estimated$condcorrs[,1]
  cor_ANZ_NAB <- tmp$stcc$Estimated$condcorrs[,2]
  cor_ANZ_WBC <- tmp$stcc$Estimated$condcorrs[,3]
  cor_CBA_NAB <- tmp$stcc$Estimated$condcorrs[,4]
  cor_CBA_WBC <- tmp$stcc$Estimated$condcorrs[,5]
  cor_NAB_WBC <- tmp$stcc$Estimated$condcorrs[,6]
  
  plot(dates,cor_ANZ_CBA,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax), axes=FALSE, lty=1,lwd=1.5)
  par(new = TRUE)
  plot(dates,cor_ANZ_NAB,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax), axes=FALSE, lty=2,lwd=1.5)
  par(new = TRUE)
  plot(dates,cor_ANZ_WBC,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax), axes=FALSE, lty=3,lwd=1.5)
  par(new = TRUE)
  plot(dates,cor_CBA_NAB,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax), axes=FALSE, lty=4,lwd=1.5)
  par(new = TRUE)
  plot(dates,cor_CBA_WBC,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax), axes=FALSE, lty=5,lwd=1.5)
  par(new = TRUE)
  plot(dates,cor_NAB_WBC,type="l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax), lty=6,lwd=1.5)
  legend("bottomright",c("ANZ-CBA","ANZ-NAB","ANZ-WBC","CBA-NAB","CBA-WBC","NAB-WBC"),lty=c(1,2,3,4,5,6),bty="n",lwd=1.5,ncol=2,y.intersp=2,inset=0.02)
  par(new = FALSE)

  dev.copy(pdf,'correlations.pdf',width=7, height=4)
  dev.off()
  
  
  
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
GARCH <- newGARCH(GARCHtype$GJR,c(0.05,0.05,0.85,0.05))
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
GARCH <- newGARCH(GARCHtype$GJR,c(0.05,0.05,0.85,0.05))
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
GARCH <- newGARCH(GARCHtype$GJR,c(0.05,0.05,0.85,0.05))
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
GARCH <- newGARCH(GARCHtype$GJR,c(0.05,0.05,0.85,0.05))
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
Tobs <- ROW(e)    

STCC <- newSTCC(Tobs)
STCC$P1 <- matrix(0.2,N,N)  # Const Corr = 0.2
diag(STCC$P1) <- 1
STCC$P2 <- matrix(0.8,N,N)  # Const Corr = 0.8
diag(STCC$P2) <- 1

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


##==============================================================##
####  Step 6: Execute multivar estimation loop  ####
##==============================================================##
####====================== Data Load =============================####
DoThis <- TRUE
if (DoThis){
  mydata <- readRDS("AUS_bankreturns_1992_2020.RDS")
  dates <- mydata$date
  e_anz <- mydata$return_anz * 100 
  e_anz <- e_anz - mean(e_anz)
  e_cba <- mydata$return_cba * 100 
  e_cba <- e_cba - mean(e_cba)
  e_nab <- mydata$return_nab * 100 
  e_nab <- e_nab - mean(e_nab)
  e_wbc <- mydata$return_wbc * 100 
  e_wbc <- e_wbc - mean(e_wbc)
}
##====================== Data Load =============================##
####====================== Object Load =============================####
if (TRUE){
  # Read AUS_4Banks manual estimation objects
  MTV1_anz <- readRDS("Output/ANZ_mtvgjr_manual.RDS")
  MTV1_cba <- readRDS("Output/CBA_mtvgjr_manual.RDS")
  MTV1_nab <- readRDS("Output/NAB_mtvgjr_manual.RDS")
  MTV1_wbc <- readRDS("Output/WBC_mtvgjr_manual.RDS")
}
####====================== Object Load =============================####

nTV <- list()
nTV[[1]] <-MTV1_anz$Estimated[[1]]$tv
nTV[[2]] <-MTV1_cba$Estimated[[1]]$tv
nTV[[3]] <-MTV1_nab$Estimated[[1]]$tv
nTV[[4]] <-MTV1_wbc$Estimated[[1]]$tv

nGARCH <- list()
nGARCH[[1]] <-MTV1_anz$Estimated[[1]]$garch
nGARCH[[2]] <-MTV1_cba$Estimated[[1]]$garch
nGARCH[[3]] <-MTV1_nab$Estimated[[1]]$garch
nGARCH[[4]] <-MTV1_wbc$Estimated[[1]]$garch


# Load up the univariate objects:
#nTV <- readRDS("AUS_4Banks_nTV.RDS")
#nGARCH <- readRDS("AUS_4Banks_nGARCH.RDS")
#STCC <- readRDS("STCC_ht1.RDS")



#Build the matrix of data:
e <- cbind(e_anz,e_cba,e_nab,e_wbc)
z_anz <- e_anz/sqrt(nTV[[1]]@g*nGARCH[[1]]@h)
z_cba <- e_cba/sqrt(nTV[[2]]@g*nGARCH[[2]]@h)
z_nab <- e_nab/sqrt(nTV[[3]]@g*nGARCH[[3]]@h)
z_wbc <- e_wbc/sqrt(nTV[[4]]@g*nGARCH[[4]]@h)
z<- cbind(z_anz,z_cba,z_nab,z_wbc)

Tobs <- NROW(e)
N <- NCOL(e)
STCC <- newSTCC(Tobs)
P1 <- matrix(0.2,N,N)
diag(P1) <- 1
STCC$P1 <- P1
P2 <- matrix(0.8,N,N)
diag(P2) <- 1
STCC$P2 <- P2
STCC$shape <- TVshape$double

stcc<-STCC
tmp <- EstimateSTCC(z,STCC)

# Set the variance targetting mode:
V_Target <- TRUE

# Set up the multivariate Return object:
mvRTN <- list()     
mvRTN$ngarch <- nGARCH
mvRTN$ntv <- nTV
mvRTN$stcc <- STCC

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
  saveRDS(m_results,"TVGJR_MultivarEstimation.RDS")

  ###  Check if things have changed!!!  ###
  
  identical(nTV,mvRTN$ntv)
  identical(nGARCH,mvRTN$ngarch)
  identical(STCC,mvRTN$stcc)
  #
}

##==============================================================##
##                            THE END
##==============================================================##



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

STCC <- newSTCC(Tobs)
STCC$P1 <- matrix(0.2,N,N)  # Const Corr = 0.2
diag(STCC$P1) <- 1
STCC$P2 <- matrix(0.8,N,N)  # Const Corr = 0.8
diag(STCC$P2) <- 1




#ntv <- mvRTN$ntv
#ngarch <- mvRTN$ngarch
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



