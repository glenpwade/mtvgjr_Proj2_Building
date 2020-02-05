try(rm(list=ls()))
gc()

##====================== Initialisation ========================####
if(T){

  #library(foreach)
  #library(doParallel)
  library(matrixcalc)
  # Multi-core math library: 
  library(RevoUtilsMath)

  #setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        #Anna's work PC - GOOGLE DRIVE
  #setwd("E:/WORK/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")   #Anna's Laptop - GOOGLE DRIVE
  setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")         # Glens Home PC
  #setwd("C:/MTVGJR_MGARCH/Project2_Building")                                              # Glens Work PC
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  ##functionsFile <- file.path(functionsPath,"functions_tvgjr_v3.r")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v7_Anna.r")
  source(functionsFile)
  
}
##====================== Initialisation ========================##


##====================== Load pre-created noise,discard & garch data & Setup DoParallel ========================##

  noiseData <- readRDS("RNorm_RefData_Sim5000_N20_T2000.RDS")
  discard_noiseData <- readRDS("RNorm_DiscardData_Sim5000_N20_T400.RDS")
  garchparsData <- readRDS("SimulationGarchPars2.RDS")
  garchparsData <- garchparsData[3:5,]  #Row 1 & 2 contain the persistance & kurtosis (not needed here)
  tvparsData <- readRDS("SimulationTVPars.RDS")
  tvparsData <- tvparsData[1,]  #Row 1 contains delta0

  # Multi-core Parallel processing:
  numcores <- 3
  Sys.setenv("MC_CORES" = numcores)
  cl <- makeCluster(numcores)
  registerDoParallel(cl, cores = numcores)
  
  
##====================== Load pre-created noise data & Setup DoParallel ========================##

  
## Manage the loops as follows:

##  1. Start with rho (re-run this file for all values of rho)
##  2. For each value of rho, write a seperate case {} for each value of N
##  3. For each N, execute a seperate 'b-loop' for each value of Tobs

##==============================================================================================##
  
  simLoops <- 5
  rho <- 0.5  # Change this and run code below:
  
# N = 2, Tobs = 25,50,100,250,500,1000,2500,5000 Garch included in data series, but not modelled
# refData <- noiseData[,1:10000]

  N <- 2
  
  # Select the first (2 x Simloops) columns from the Noise Data to use here:
  refData <- noiseData[,1:1000]
  discardData <- discard_noiseData[,1:100]
  
  #saveAs <- paste0("Simulated_Distributions/ccDist_equi0",100*rho,"_N",N,"_IgnoreGarch2.1.RDS")
  ccDistResults <- NULL
  
  # CCC
  CCC <- list()
  CCC$P <- matrix(rho,nrow=N,ncol=N)
  diag(CCC$P) <- 1
  
  CORR <- list()
  CORR$type <- "CCC"
  CORR$CCC <- CCC
  
  Tobs <- 250
  Testing_TV <- TRUE
  
  nGarchData <- list()
  nTVData <- list()
  
if(Testing_TV){

  for (n in 1:10000) {
    GARCH <- list()
    GARCH$type <- GARCHtype$none
    GARCH$pars <- c(0.05,0.05,0.95)
    nGarchData[[n]] <- GARCH
    #
    TV <- list()
    TV$var_target <- TRUE
    TV$Tobs <- Tobs
    TV$shape[1] <- TRshape$none
    tvparsData[n] <- 0
    TV$delta0 <- tvparsData[n]
    TV$Estimated$delta0 <- tvparsData[n]
    TV$Estimated$lastdelta0 <- tvparsData[n]
    TV$condvars <- tvparsData[n]
    TV$Estimated$condvars <- tvparsData[n]
    nTVData[[n]] <- TV
  }
} else
  {

    for (n in 1:10000) {
      GARCH <- list()
      GARCH$type <- GARCHtype$constant
      garchparsData[,n] <- c(0.05,0.05,0.95)
      GARCH$pars <- garchparsData[,n]
      nGarchData[[n]] <- GARCH
      #
      TV <- list()
      TV$var_target <- FALSE
      TV$Tobs <- Tobs
      TV$shape[1] <- TRshape$none
      TV$delta0 <- 0
      TV$Estimated$delta0 <- 0
      TV$Estimated$lastdelta0 <- 0
      nTVData[[n]] <- TV
    }
  }
  
  if(T) {

    # STCC - used for H1
    STCC <- list()
    STCC$st <- (1:Tobs)/Tobs
    # Simulate the distribution:
    
      b <- 1
      end <- b*N
      start <- (end - N) + 1
      refData1 <- refData[(1:Tobs),(start:end)]
      discardData1 <- discardData[(1:10),(start:end)]

      if(Testing_TV){
        # Generate Simulation Data with specific CORR & GARCH
        nTV <- nTVData[start:end]
        e <- GenerateRefData(corr=CORR,ngarch=NULL,ntv=nTV,e=refData1,e_discard=discardData1,numseries=N,tobs=Tobs,saveas=NULL)  
      }else{
        ## Generate Simulation Data with specific CORR & TV
        nGarch <- nGarchData[start:end]
        e <- GenerateRefData(corr=CORR,ngarch=nGarch,ntv=NULL,e=refData1,e_discard=discardData1,numseries=N,tobs=Tobs,saveas=NULL)
        
      }
      
      z <- e
      
      if(Testing_TV){
         for (n in 1:N) {
        #   TV <- list()
        #   GARCH <- list()
        #   TV$var_target <- TRUE
        #   TV$shape <- TRshape$none
        #   TV$Estimated$delta0 <- var(e[,n])
        #   TV$Estimated$condvars <- rep(var(e[,n]),Tobs)
           nTV[[n]]$condvars <- rep(var(e[,n]),Tobs)
           z[,n] <- e[,n]/sqrt(var(e[,n]))
           nTV[[n]]$Estimated$condvars <- nTV[[n]]$condvars
        #   nTV[[n]] <- TV
           nGarch[[n]]$type <- GARCHtype$none
        #   GARCH$condvars <- rep(1,Tobs)
        #   nGarch[[n]] <- GARCH
         }
        nTV$Type <- TRshape$single
        nGarch$Type <- GARCHtype$none
      } else {
         for (n in 1:N) {
        #   TV <- list()
        #   GARCH <- list()
        #   GARCH$type <- GARCHtype$constant
        #   GARCH$pars <- c(0.02,0.03,0.95)
           nGarch[[n]]$condvars <- rep(var(e[,n]),Tobs)
           z[,n] <- e[,n]/sqrt(var(e[,n]))
        #   nGarch[[n]] <- GARCH
        #   TV$delta0 <- 0
        #   TV$var_target <- TRUE
        #   TV$shape <- TRshape$none
        #   TV$condvars <- 0
        #   nTV[[n]] <- TV
         }
        nTV$Type <- TRshape$none
        nGarch$Type <- GARCHtype$constant
      }
      
      CCC <- list()
      CCC$P <- cor(z)

      H0 <- list()
      H0$CCC<- CCC
      H0$nGARCH <- nGarch
      H0$nTV <- nTV
      H1 <- STCC

      T1 <- myTest.CCCvSTCC.LM(e=e,H0=H0,H1=H1,testorder = 1)
      T2 <- myTest.CCCvSTCC.LM(e=e,H0=H0,H1=H1,testorder = 2)
      #Return
      c(T1,T2)

    
    ccDistResults <- cbind(ccDistResults,ccDist)

  } #End: T=25
  
## ========================  Close DoParallel & clean up ========================  ##
  try(unregisterDoParallel(cl))
  try(rm(list=ls()))
  try(gc())
  
  
##====================== End: Simulate Probability Distribution ===========================##
stop()
  
  
##====================== Utility Section ===========================##
  
  ccDistResults <- readRDS("ccDist_equi025_N3_IgnoreGarch2.RDS")
  
  # Calculate the 1,5 & 10% pvalue sizes:
  tmp = 0
  sum(ccDistResults[,tmp]>qchisq(p=0.01,df=N*(N-1)/2,lower.tail=FALSE))/length(ccDistResults[,tmp])
  sum(ccDistResults[,tmp]>qchisq(p=0.05,df=N*(N-1)/2,lower.tail=FALSE))/length(ccDistResults[,tmp])
  sum(ccDistResults[,tmp]>qchisq(p=0.10,df=N*(N-1)/2,lower.tail=FALSE))/length(ccDistResults[,tmp])
  sum(ccDistResults[,(tmp+1)]>qchisq(p=0.01,df=2*N*(N-1)/2,lower.tail=FALSE))/length(ccDistResults[,(tmp+1)])
  sum(ccDistResults[,(tmp+1)]>qchisq(p=0.05,df=2*N*(N-1)/2,lower.tail=FALSE))/length(ccDistResults[,(tmp+1)])
  sum(ccDistResults[,(tmp+1)]>qchisq(p=0.10,df=2*N*(N-1)/2,lower.tail=FALSE))/length(ccDistResults[,(tmp+1)])

##====================== Utility Section ===========================##


