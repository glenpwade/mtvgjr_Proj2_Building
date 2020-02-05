##====================== Initialisation ========================##

# We want to calculate the optimum Block length for the Block-Bootstrap

rm(list=ls())
gc(T)

if (TRUE){
  library(graphics)
  library(foreach)
  library(doParallel)
  library(moments)
  library(Matrix)
  
  #setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's work PC - GOOGLE DRIVE
  #setwd("~/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Anna's laptop - GOOGLE DRIVE
  setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")        # Glen's work PC - GOOGLE DRIVE
  
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v6.r")
  source(functionsFile)
  
  ## Setup the parallel backend environment ##
  numcores <- 3
  Sys.setenv("MC_CORES" = numcores)
  cl <- makeCluster(numcores)
  registerDoParallel(cl, cores = numcores)
  
  ## Remove parallel backend - if already setup ##
  if(FALSE) {
    unregisterDoParallel(cl)
    rm(cl)  
  }
  
  
}


##====================== Initialisation ========================##


#### ====================== Data Setup ============================ ####
DoThis <- TRUE
if (DoThis){
  
G1 <- newGARCH(GARCHtype$general,c(0.03,0.07,0.9))
simData <- GenerateRefData_WithGarch(garch=G1,tsamples = 1000)

}
##====================== Data Setup ============================##


####====================== BootStrap =============================####

N <- 10
numBS <- 99
testStatDist <- matrix(NA,nrow=numBS+1,ncol=N)
set.seed(1999)

## Outer Loop: Want to return 1000 column Vectors, where 1st element is the RefTest stat & next 999 are BS
testStatDist <- foreach(b = 1:N, .inorder=FALSE, .combine=cbind, .verbose = FALSE) %dopar% {

BS_Result <- matrix(nrow = 1,ncol = numBS+1)

#----  Get a Reference Statistic:  ----#
e <- simData[,1]
# Create a TV object with just delta0 & st:
ST <- 1:length(e)/length(e)
TV <- newTV(del0 = 1,st = ST)
# Estimate the univar series:
TV <- EstimateTV(e,TV)
# Quick look at the Estimation Outputs:
if (FALSE) {
  TV$Estimated$value
  TV$Estimated$delta0  
  plot(e,type="l")
}

# RefTest Stat in 1st position
BS_Result[1,1] <- myTest.TV.noGARCH.TR2(e,TV)

## Hack: Take this out of the loop below - just for performance
# Standardise TV
ST <- 1:990/990
TV$st <- ST
TV$Tobs <- 990
TV$Estimated$condvars <- 1.0

## ----  Now the inner loop  ---- ##
for (b in 2:numBS){
  set.seed(b*1009)
  
  e <- simData[,b]
  B_len <- 90  # Calculated based on T=1000
  # T <- length(e)
  # numBlocks <- round(T/B_len)
  # T <- numBlocks * B_len
  
  #  Hack for performance:
  T <- 990
  numBlocks <- 11
  
  randDraw <- round(runif(numBlocks,min=0,max=numBlocks-1))
  e_b <- vector("numeric",0)
  
  for (i in seq_along(randDraw)) {
    e_b <- c(e_b,e[(randDraw[i]*B_len+1):(randDraw[i]*B_len+B_len)])
  }
  
  # # Standardise TV
  # ST <- 1:NROW(e_b)/NROW(e_b)
  # TV$st <- ST
  # TV$Tobs <- NROW(e_b)
  # TV$Estimated$condvars <- 1.0
  #Return:
  try(BS_Result[1,b] <- myTest.TV.noGARCH.TR2(e_b,TV))
  
  
} # End: for(b in 2:numBS,...
BS_Result

}

# Plot the Test Stat distribution:
hist(testStatDist[1,],30)

# Calculate p-value:
pVal <- length(testStatDist[1,which(refTest > testStatDist[,1])])/(numBS+1)


##====================== Bootstrap ============================##

##==============================================================##
##                            THE END
##==============================================================##


