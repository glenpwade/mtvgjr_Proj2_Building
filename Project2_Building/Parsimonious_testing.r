try(rm(list=ls()))

##====================== Initialisation ========================##
if(T){
  
  library(foreach)
  library(matrixcalc)
  library(RevoUtilsMath)
  RevoUtilsMath::setMKLthreads(3)
  
  
  setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project10_Eigentest")        #Anna's work PC - GOOGLE DRIVE
  #setwd("E:/WORK/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project10_Eigentest")   #Anna's Laptop - GOOGLE DRIVE
  #setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")    # Glen's home-office PC
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v9.r")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"clsMTVGARCH.R")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"Test_CCCParsim_LM.R")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"Test_CCCvSTCC_LM.R")
  source(functionsFile)
  
}


##====================== Load pre-created noise,discard & garch data & Setup DoParallel ========================##
if(T){
  #noiseData <- readRDS("RNorm_RefData_Sim5000_N20_T2000.RDS")
  #discard_noiseData <- readRDS("RNorm_DiscardData_Sim5000_N20_T400.RDS")
  noiseData <- readRDS("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building/RNorm_RefData_Sim5000_N20_T5000.RDS")
  discard_noiseData <- readRDS("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building/RNorm_DiscardData_Sim5000_N20_T1500.RDS")
}

##====================== Load pre-created noise data & Setup DoParallel ========================##

## ===============  START ================ ####  
ccDistResults <- NULL

# Change the values below and run code below again:
rho <- 1/3
GARCH_InRefData <- TRUE
per <- 0.95
kur <- 4
GARCH_InTest <- TRUE
N <- 2
Tobs <- 1000
Bobs <- 1000
Bobs_start <- 1
Bobs_end <- Bobs
refData <- noiseData[(1:Tobs),(1:(N*Bobs))]
discardData <- discard_noiseData[,(1:(N*Bobs))]
#fName <- paste0("ccDist_CEC",round(rho*100),"_N",N,"_T",Tobs,"_ConstTV_GARCH_p",per*100,"_k",kur,".RDS")

cat("\nStarting Loop: N =",N, "  #Obs:",Tobs)

CCC <- list()
#CCC$P <- matrix(rho,nrow=N,ncol=N)
## CEC ##
if(F){
  CCC$P <- matrix(rho,nrow=N,ncol=N)
  diag(CCC$P) <- 1
}
## CCC Toepliz ##
if(T){
  CCC$P <- toeplitz(rho^seq.int(0, N-1))
}

CORR <- list()
CORR$type <- "CCC"
CORR$CCC <- CCC

nTVgenerate <- list()
nGarchgenerate <- list()
nTV <- list()
nGarch <- list()

# Create TV & GARCH objects to be used in Generating RefData  
# Create TV & GARCH objects to be used in Generating RefData  
st <- (1:Tobs)/Tobs
TV <- tv(st,tvshape$delta0only)

for (n in 1:N) {
  if(GARCH_InRefData){
    MTVG1 <- mtvgarch(TV,garchtype$general)
    alphapar <- sqrt((1-per^2)*(1-3/kur)*0.5)
    betapar <- per-alphapar
    omegapar <- 1-alphapar-betapar
    MTVG1$garch$pars['omega',1] <- omegapar
    MTVG1$garch$pars['alpha',1] <- alphapar
    MTVG1$garch$pars['beta',1] <- betapar
  } else {
    MTVG1 <- mtvgarch(TV,garchtype$noGarch)
    MTVG1$tv$delta0 <- n+1
  }
  nTVgenerate[[n]] <- MTVG1$tv
  nGarchgenerate[[n]] <- MTVG1$garch
}

# STCC - used for H1
STCC <- list()
STCC$st <- (1:Tobs)/Tobs
# Simulate the distribution:
b=1

for (b in Bobs_start:Bobs_end){
  #b=5
  #ccDist <- foreach(b = 1:simLoops, .inorder=FALSE, .combine=rbind, .verbose = FALSE) %do% {
  
  end <- b*N
  start <- (end - N) + 1
  refData1 <- refData[(1:Tobs),(start:end)]
  discardData1 <- discardData[(1:50),(start:end)]
  
  # Generate Simulation Data with specific CORR & GARCH
  e <- GenerateRefData(corr=CORR,ngarch=nGarchgenerate,ntv=nTVgenerate,e=refData1,e_discard=discardData1,numseries=N,tobs=Tobs,saveas=NULL)  
  
  w <- e
  z <- w
  
  #plot(e[,1],type='l')
  #plot(e[,2],type='l')
  
  # Create TV & GARCH objects to be used in Estimation and Test  
  for (n in 1:N) {
    st <- (1:Tobs)/Tobs
    TV <- tv(st,tvshape$delta0only)
    TV <- estimateTV(e[,n],TV)
    if(GARCH_InTest){
      MTVG2 <- mtvgarch(TV,garchtype$general)
      MTVG2$garch$pars['omega',1] <- 0.03
      MTVG2$garch$pars['alpha',1] <- 0.02
      MTVG2$garch$pars['beta',1] <- 0.95
    } else {
      MTVG2 <- mtvgarch(TV,garchtype$noGarch)
    }
    nTV[[n]] <- MTVG2$tv
    nGarch[[n]] <- MTVG2$garch
  }  
  
  for (n in 1:N) {
    w[,n] <- e[,n]/sqrt(nTV[[n]]@g)
    if(nGarch[[n]]$type != garchtype$noGarch){
      nGarch[[n]]$optimcontrol$reltol <- 1e-4
      nGarch[[n]]$optimcontrol$parscale <- c(1,1,3)
      nGarch[[n]]$optimcontrol$ndeps <- c(1e-3,1e-3,1e-3)
      nGarch[[n]] <- estimateGARCH(w[,n],nGarch[[n]])  
      z[,n] <- w[,n]/sqrt(nGarch[[n]]@h)
    }
  }
  
  #plot(nTV[[n]]@g,typ='l')
  var(z)
  
  CCC <- list()
  CCC$P <- cor(z)
  
  H0 <- list()
  H0$CCC <- CCC
  H0$nGARCH <- nGarch
  H0$nTV <- nTV
  H1 <- STCC
  
  T1<-T2<-T3<-T4<-T5<-0
  
  T1 <- myTest.CCCvSTCC.LM(e=e,H0=H0,H1=H1,testorder = 1)
  T2 <- myTest.CCCvSTCC.LM(e=e,H0=H0,H1=H1,testorder = 2)
  T3 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 1)
  T4 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 2)
  #T5 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 3)
  #Return
  c(N,Tobs,T1,pchisq(T1, df=N*(N-1)/2 , lower.tail=FALSE),T2,pchisq(T2, df=2*N*(N-1)/2 , lower.tail=FALSE),T3,pchisq(T3, df=N-1 , lower.tail=FALSE),T4,pchisq(T4, df=2*(N-1), lower.tail=FALSE))
  ccDist <- c(b,N,Tobs,T1,T2,T3,T4)
  # ccDist <- c(b,N,Tobs,T1,T2,T3,T4,T5)
  
  ccDistResults <- rbind(ccDistResults,ccDist)
  # df1 <- N*(N-1)/2
  # df2 <- 2*N*(N-1)/2
  # df3 <- N-1
  # df4 <- 2*(N-1)
  # df5 <- 3*(N-1)
  # print(c(b,":T1",round(sum(ccDistResults[,4]>qchisq(.95,df=df1))/b,4),"T2",round(sum(ccDistResults[,5]>qchisq(.95,df=df2))/b,4),"T3",round(sum(ccDistResults[,6]>qchisq(.95,df=df3))/b,4),"T4",round(sum(ccDistResults[,7]>qchisq(.95,df=df4))/b,4)),quote=FALSE)
}  # End of b 


colnames(ccDistResults) <- c("b","N","Tobs","T1","T2","T3","T4")
# Chi-Squared(df=testorder*N(N-1)/2)
sum(ccDistResults[,"T1"]>qchisq(.99, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T1"]>qchisq(.95, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T1"]>qchisq(.90, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T2"]>qchisq(.99, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T2"]>qchisq(.95, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T2"]>qchisq(.90, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T3"]>qchisq(.99, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T3"]>qchisq(.95, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T3"]>qchisq(.90, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T4"]>qchisq(.99, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
sum(ccDistResults[,"T4"]>qchisq(.95, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
sum(ccDistResults[,"T4"]>qchisq(.90, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
# sum(ccDistResults[,"T5"]>qchisq(.99, df=3*(N-1)))/NROW(ccDistResults[,"T5"])
# sum(ccDistResults[,"T5"]>qchisq(.95, df=3*(N-1)))/NROW(ccDistResults[,"T5"])
# sum(ccDistResults[,"T5"]>qchisq(.90, df=3*(N-1)))/NROW(ccDistResults[,"T5"])

#saveRDS(ccDistResults,fName)









try(rm(list=ls()))

##====================== Initialisation ========================##
if(T){
  
  library(foreach)
  library(matrixcalc)
  library(RevoUtilsMath)
  RevoUtilsMath::setMKLthreads(3)
  
  
  setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project10_Eigentest")        #Anna's work PC - GOOGLE DRIVE
  #setwd("E:/WORK/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project10_Eigentest")   #Anna's Laptop - GOOGLE DRIVE
  #setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")    # Glen's home-office PC
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v9.r")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"clsMTVGARCH.R")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"Test_CCCParsim_LM.R")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"Test_CCCvSTCC_LM.R")
  source(functionsFile)
  
}


##====================== Load pre-created noise,discard & garch data & Setup DoParallel ========================##
if(T){
  #noiseData <- readRDS("RNorm_RefData_Sim5000_N20_T2000.RDS")
  #discard_noiseData <- readRDS("RNorm_DiscardData_Sim5000_N20_T400.RDS")
  noiseData <- readRDS("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building/RNorm_RefData_Sim5000_N20_T5000.RDS")
  discard_noiseData <- readRDS("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building/RNorm_DiscardData_Sim5000_N20_T1500.RDS")
}

##====================== Load pre-created noise data & Setup DoParallel ========================##

## ===============  START ================ ####  
ccDistResults <- NULL

# Change the values below and run code below again:
rho <- 1/3
GARCH_InRefData <- TRUE
per <- 0.95
kur <- 4
GARCH_InTest <- TRUE
N <- 50
Tobs <- 1000
Bobs <- 1000
Bobs_start <- 1
Bobs_end <- Bobs
refData <- noiseData[(1:Tobs),(1:(N*Bobs))]
discardData <- discard_noiseData[,(1:(N*Bobs))]
#fName <- paste0("ccDist_CEC",round(rho*100),"_N",N,"_T",Tobs,"_ConstTV_GARCH_p",per*100,"_k",kur,".RDS")

cat("\nStarting Loop: N =",N, "  #Obs:",Tobs)

CCC <- list()
#CCC$P <- matrix(rho,nrow=N,ncol=N)
## CEC ##
if(F){
  CCC$P <- matrix(rho,nrow=N,ncol=N)
  diag(CCC$P) <- 1
}
## CCC Toepliz ##
if(T){
  CCC$P <- toeplitz(rho^seq.int(0, N-1))
}

CORR <- list()
CORR$type <- "CCC"
CORR$CCC <- CCC

nTVgenerate <- list()
nGarchgenerate <- list()
nTV <- list()
nGarch <- list()

# Create TV & GARCH objects to be used in Generating RefData  
# Create TV & GARCH objects to be used in Generating RefData  
st <- (1:Tobs)/Tobs
TV <- tv(st,tvshape$delta0only)

for (n in 1:N) {
  if(GARCH_InRefData){
    MTVG1 <- mtvgarch(TV,garchtype$general)
    alphapar <- sqrt((1-per^2)*(1-3/kur)*0.5)
    betapar <- per-alphapar
    omegapar <- 1-alphapar-betapar
    MTVG1$garch$pars['omega',1] <- omegapar
    MTVG1$garch$pars['alpha',1] <- alphapar
    MTVG1$garch$pars['beta',1] <- betapar
  } else {
    MTVG1 <- mtvgarch(TV,garchtype$noGarch)
    MTVG1$tv$delta0 <- n+1
  }
  nTVgenerate[[n]] <- MTVG1$tv
  nGarchgenerate[[n]] <- MTVG1$garch
}

# STCC - used for H1
STCC <- list()
STCC$st <- (1:Tobs)/Tobs
# Simulate the distribution:
b=1

for (b in Bobs_start:Bobs_end){
  #b=5
  #ccDist <- foreach(b = 1:simLoops, .inorder=FALSE, .combine=rbind, .verbose = FALSE) %do% {
  
  end <- b*N
  start <- (end - N) + 1
  refData1 <- refData[(1:Tobs),(start:end)]
  discardData1 <- discardData[(1:50),(start:end)]
  
  # Generate Simulation Data with specific CORR & GARCH
  e <- GenerateRefData(corr=CORR,ngarch=nGarchgenerate,ntv=nTVgenerate,e=refData1,e_discard=discardData1,numseries=N,tobs=Tobs,saveas=NULL)  
  
  w <- e
  z <- w
  
  #plot(e[,1],type='l')
  #plot(e[,2],type='l')
  
  # Create TV & GARCH objects to be used in Estimation and Test  
  for (n in 1:N) {
    st <- (1:Tobs)/Tobs
    TV <- tv(st,tvshape$delta0only)
    TV <- estimateTV(e[,n],TV)
    if(GARCH_InTest){
      MTVG2 <- mtvgarch(TV,garchtype$general)
      MTVG2$garch$pars['omega',1] <- 0.03
      MTVG2$garch$pars['alpha',1] <- 0.02
      MTVG2$garch$pars['beta',1] <- 0.95
    } else {
      MTVG2 <- mtvgarch(TV,garchtype$noGarch)
    }
    nTV[[n]] <- MTVG2$tv
    nGarch[[n]] <- MTVG2$garch
  }  
  
  for (n in 1:N) {
    w[,n] <- e[,n]/sqrt(nTV[[n]]@g)
    if(nGarch[[n]]$type != garchtype$noGarch){
      nGarch[[n]]$optimcontrol$reltol <- 1e-4
      nGarch[[n]]$optimcontrol$parscale <- c(1,1,3)
      nGarch[[n]]$optimcontrol$ndeps <- c(1e-3,1e-3,1e-3)
      nGarch[[n]] <- estimateGARCH(w[,n],nGarch[[n]])  
      z[,n] <- w[,n]/sqrt(nGarch[[n]]@h)
    }
  }
  
  #plot(nTV[[n]]@g,typ='l')
  var(z)
  
  CCC <- list()
  CCC$P <- cor(z)
  
  H0 <- list()
  H0$CCC <- CCC
  H0$nGARCH <- nGarch
  H0$nTV <- nTV
  H1 <- STCC
  
  T1<-T2<-T3<-T4<-T5<-0
  
  T1 <- myTest.CCCvSTCC.LM(e=e,H0=H0,H1=H1,testorder = 1)
  T2 <- myTest.CCCvSTCC.LM(e=e,H0=H0,H1=H1,testorder = 2)
  T3 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 1)
  T4 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 2)
  #T5 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 3)
  #Return
  c(N,Tobs,T1,pchisq(T1, df=N*(N-1)/2 , lower.tail=FALSE),T2,pchisq(T2, df=2*N*(N-1)/2 , lower.tail=FALSE),T3,pchisq(T3, df=N-1 , lower.tail=FALSE),T4,pchisq(T4, df=2*(N-1), lower.tail=FALSE))
  ccDist <- c(b,N,Tobs,T1,T2,T3,T4)
  # ccDist <- c(b,N,Tobs,T1,T2,T3,T4,T5)
  
  ccDistResults <- rbind(ccDistResults,ccDist)
  df1 <- N*(N-1)/2
  df2 <- 2*N*(N-1)/2
  df3 <- N-1
  df4 <- 2*(N-1)
  df5 <- 3*(N-1)
  print(c(b,":T1",round(sum(ccDistResults[,4]>qchisq(.95,df=df1))/b,4),"T2",round(sum(ccDistResults[,5]>qchisq(.95,df=df2))/b,4),"T3",round(sum(ccDistResults[,6]>qchisq(.95,df=df3))/b,4),"T4",round(sum(ccDistResults[,7]>qchisq(.95,df=df4))/b,4),"T5",round(sum(ccDistResults[,8]>qchisq(.95,df=df5))/b,4)),quote=FALSE)
}  # End of b 


colnames(ccDistResults) <- c("b","N","Tobs","T1","T2","T3","T4","T5")
# Chi-Squared(df=testorder*N(N-1)/2)
sum(ccDistResults[,"T1"]>qchisq(.99, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T1"]>qchisq(.95, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T1"]>qchisq(.90, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T2"]>qchisq(.99, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T2"]>qchisq(.95, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T2"]>qchisq(.90, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T3"]>qchisq(.99, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T3"]>qchisq(.95, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T3"]>qchisq(.90, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T4"]>qchisq(.99, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
sum(ccDistResults[,"T4"]>qchisq(.95, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
sum(ccDistResults[,"T4"]>qchisq(.90, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
# sum(ccDistResults[,"T5"]>qchisq(.99, df=3*(N-1)))/NROW(ccDistResults[,"T5"])
# sum(ccDistResults[,"T5"]>qchisq(.95, df=3*(N-1)))/NROW(ccDistResults[,"T5"])
# sum(ccDistResults[,"T5"]>qchisq(.90, df=3*(N-1)))/NROW(ccDistResults[,"T5"])

#saveRDS(ccDistResults,fName)






try(rm(list=ls()))

##====================== Initialisation ========================##
if(T){
  
  library(foreach)
  library(matrixcalc)
  library(RevoUtilsMath)
  RevoUtilsMath::setMKLthreads(3)
  
  
  setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project10_Eigentest")        #Anna's work PC - GOOGLE DRIVE
  #setwd("E:/WORK/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project10_Eigentest")   #Anna's Laptop - GOOGLE DRIVE
  #setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")    # Glen's home-office PC
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v9.r")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"clsMTVGARCH.R")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"Test_CCCParsim_LM.R")
  source(functionsFile)
  functionsFile <- file.path(functionsPath,"Test_CCCvSTCC_LM.R")
  source(functionsFile)
  
}


##====================== Load pre-created noise,discard & garch data & Setup DoParallel ========================##
if(T){
  #noiseData <- readRDS("RNorm_RefData_Sim5000_N20_T2000.RDS")
  #discard_noiseData <- readRDS("RNorm_DiscardData_Sim5000_N20_T400.RDS")
  noiseData <- readRDS("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building/RNorm_RefData_Sim5000_N20_T5000.RDS")
  discard_noiseData <- readRDS("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building/RNorm_DiscardData_Sim5000_N20_T1500.RDS")
}

##====================== Load pre-created noise data & Setup DoParallel ========================##

## ===============  START ================ ####  
ccDistResults <- NULL

# Change the values below and run code below again:
rho <- 1/3
GARCH_InRefData <- TRUE
per <- 0.95
kur <- 4
GARCH_InTest <- TRUE
N <- 100
Tobs <- 1000
Bobs <- 1000
Bobs_start <- 1
Bobs_end <- Bobs
refData <- noiseData[(1:Tobs),(1:(N*Bobs))]
discardData <- discard_noiseData[,(1:(N*Bobs))]
#fName <- paste0("ccDist_CEC",round(rho*100),"_N",N,"_T",Tobs,"_ConstTV_GARCH_p",per*100,"_k",kur,".RDS")

cat("\nStarting Loop: N =",N, "  #Obs:",Tobs)

CCC <- list()
#CCC$P <- matrix(rho,nrow=N,ncol=N)
## CEC ##
if(F){
  CCC$P <- matrix(rho,nrow=N,ncol=N)
  diag(CCC$P) <- 1
}
## CCC Toepliz ##
if(T){
  CCC$P <- toeplitz(rho^seq.int(0, N-1))
}

CORR <- list()
CORR$type <- "CCC"
CORR$CCC <- CCC

nTVgenerate <- list()
nGarchgenerate <- list()
nTV <- list()
nGarch <- list()

# Create TV & GARCH objects to be used in Generating RefData  
# Create TV & GARCH objects to be used in Generating RefData  
st <- (1:Tobs)/Tobs
TV <- tv(st,tvshape$delta0only)

for (n in 1:N) {
  if(GARCH_InRefData){
    MTVG1 <- mtvgarch(TV,garchtype$general)
    alphapar <- sqrt((1-per^2)*(1-3/kur)*0.5)
    betapar <- per-alphapar
    omegapar <- 1-alphapar-betapar
    MTVG1$garch$pars['omega',1] <- omegapar
    MTVG1$garch$pars['alpha',1] <- alphapar
    MTVG1$garch$pars['beta',1] <- betapar
  } else {
    MTVG1 <- mtvgarch(TV,garchtype$noGarch)
    MTVG1$tv$delta0 <- n+1
  }
  nTVgenerate[[n]] <- MTVG1$tv
  nGarchgenerate[[n]] <- MTVG1$garch
}

# STCC - used for H1
STCC <- list()
STCC$st <- (1:Tobs)/Tobs
# Simulate the distribution:
b=1

for (b in Bobs_start:Bobs_end){
  #b=5
  #ccDist <- foreach(b = 1:simLoops, .inorder=FALSE, .combine=rbind, .verbose = FALSE) %do% {
  
  end <- b*N
  start <- (end - N) + 1
  refData1 <- refData[(1:Tobs),(start:end)]
  discardData1 <- discardData[(1:50),(start:end)]
  
  # Generate Simulation Data with specific CORR & GARCH
  e <- GenerateRefData(corr=CORR,ngarch=nGarchgenerate,ntv=nTVgenerate,e=refData1,e_discard=discardData1,numseries=N,tobs=Tobs,saveas=NULL)  
  
  w <- e
  z <- w
  
  #plot(e[,1],type='l')
  #plot(e[,2],type='l')
  
  # Create TV & GARCH objects to be used in Estimation and Test  
  for (n in 1:N) {
    st <- (1:Tobs)/Tobs
    TV <- tv(st,tvshape$delta0only)
    TV <- estimateTV(e[,n],TV)
    if(GARCH_InTest){
      MTVG2 <- mtvgarch(TV,garchtype$general)
      MTVG2$garch$pars['omega',1] <- 0.03
      MTVG2$garch$pars['alpha',1] <- 0.02
      MTVG2$garch$pars['beta',1] <- 0.95
    } else {
      MTVG2 <- mtvgarch(TV,garchtype$noGarch)
    }
    nTV[[n]] <- MTVG2$tv
    nGarch[[n]] <- MTVG2$garch
  }  
  
  for (n in 1:N) {
    w[,n] <- e[,n]/sqrt(nTV[[n]]@g)
    if(nGarch[[n]]$type != garchtype$noGarch){
      nGarch[[n]]$optimcontrol$reltol <- 1e-4
      nGarch[[n]]$optimcontrol$parscale <- c(1,1,3)
      nGarch[[n]]$optimcontrol$ndeps <- c(1e-3,1e-3,1e-3)
      nGarch[[n]] <- estimateGARCH(w[,n],nGarch[[n]])  
      z[,n] <- w[,n]/sqrt(nGarch[[n]]@h)
    }
  }
  
  #plot(nTV[[n]]@g,typ='l')
  var(z)
  
  CCC <- list()
  CCC$P <- cor(z)
  
  H0 <- list()
  H0$CCC <- CCC
  H0$nGARCH <- nGarch
  H0$nTV <- nTV
  H1 <- STCC
  
  T1<-T2<-T3<-T4<-T5<-0
  
  T1 <- myTest.CCCvSTCC.LM(e=e,H0=H0,H1=H1,testorder = 1)
  T2 <- myTest.CCCvSTCC.LM(e=e,H0=H0,H1=H1,testorder = 2)
  T3 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 1)
  T4 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 2)
  #T5 <- myTest.CCCParsim.LM.new(e=e,H0=H0,H1=H1,testorder = 3)
  #Return
  c(N,Tobs,T1,pchisq(T1, df=N*(N-1)/2 , lower.tail=FALSE),T2,pchisq(T2, df=2*N*(N-1)/2 , lower.tail=FALSE),T3,pchisq(T3, df=N-1 , lower.tail=FALSE),T4,pchisq(T4, df=2*(N-1), lower.tail=FALSE))
  ccDist <- c(b,N,Tobs,T1,T2,T3,T4)
  # ccDist <- c(b,N,Tobs,T1,T2,T3,T4,T5)
  
  ccDistResults <- rbind(ccDistResults,ccDist)
  df1 <- N*(N-1)/2
  df2 <- 2*N*(N-1)/2
  df3 <- N-1
  df4 <- 2*(N-1)
  df5 <- 3*(N-1)
  print(c(b,":T1",round(sum(ccDistResults[,4]>qchisq(.95,df=df1))/b,4),"T2",round(sum(ccDistResults[,5]>qchisq(.95,df=df2))/b,4),"T3",round(sum(ccDistResults[,6]>qchisq(.95,df=df3))/b,4),"T4",round(sum(ccDistResults[,7]>qchisq(.95,df=df4))/b,4),"T5",round(sum(ccDistResults[,8]>qchisq(.95,df=df5))/b,4)),quote=FALSE)
}  # End of b 


colnames(ccDistResults) <- c("b","N","Tobs","T1","T2","T3","T4","T5")
# Chi-Squared(df=testorder*N(N-1)/2)
sum(ccDistResults[,"T1"]>qchisq(.99, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T1"]>qchisq(.95, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T1"]>qchisq(.90, df=1*N*(N-1)/2))/NROW(ccDistResults[,"T1"])
sum(ccDistResults[,"T2"]>qchisq(.99, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T2"]>qchisq(.95, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T2"]>qchisq(.90, df=2*N*(N-1)/2))/NROW(ccDistResults[,"T2"])
sum(ccDistResults[,"T3"]>qchisq(.99, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T3"]>qchisq(.95, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T3"]>qchisq(.90, df=1*(N-1)))/NROW(ccDistResults[,"T3"])
sum(ccDistResults[,"T4"]>qchisq(.99, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
sum(ccDistResults[,"T4"]>qchisq(.95, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
sum(ccDistResults[,"T4"]>qchisq(.90, df=2*(N-1)))/NROW(ccDistResults[,"T4"])
# sum(ccDistResults[,"T5"]>qchisq(.99, df=3*(N-1)))/NROW(ccDistResults[,"T5"])
# sum(ccDistResults[,"T5"]>qchisq(.95, df=3*(N-1)))/NROW(ccDistResults[,"T5"])
# sum(ccDistResults[,"T5"]>qchisq(.90, df=3*(N-1)))/NROW(ccDistResults[,"T5"])

#saveRDS(ccDistResults,fName)



