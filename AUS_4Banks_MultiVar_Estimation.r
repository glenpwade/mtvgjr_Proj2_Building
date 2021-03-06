
##  ============ Initialisation ==========  ####
rm(list=ls())
gc(T)

if (T){
  library(graphics)
  library(foreach)
  library(doParallel)
  RevoUtilsMath::setMKLthreads(4)

  #setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building_2020")  # Anna's work PC - GOOGLE DRIVE
  #setwd("D:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building_2020")   # Anna's laptop - GOOGLE DRIVE
  setwd("C:/Source/MTVGJR")
  
  source("clsCORR.r")
  source("clsTV.r")
  source("clsGARCH.r")

}
##====================== Initialisation ========================##

##====================== Data Load ============================####
if (T){
  # Read AUS_4Banks data from saved file
  mydata <- readRDS("AUS_bankreturns_1992_2020.RDS")
  dates <- mydata$date
  
  # Create scaled-up & de-meaned individual data series vectors (Percentage Returns):
  e_anz <- mydata$return_anz * 100 
  e_anz <- e_anz - mean(e_anz)
  
  e_cba <- mydata$return_cba * 100 
  e_cba <- e_cba - mean(e_cba)
  
  e_nab <- mydata$return_nab * 100 
  e_nab <- e_nab - mean(e_nab)
  
  e_wbc <- mydata$return_wbc * 100 
  e_wbc <- e_wbc - mean(e_wbc)
  
  MTV_anz <- readRDS("Output/ANZ_mtvgjr_manual.RDS")
  MTV_cba <- readRDS("Output/CBA_mtvgjr_manual.RDS")
  MTV_nab <- readRDS("Output/NAB_mtvgjr_manual.RDS")
  MTV_wbc <- readRDS("Output/WBC_mtvgjr_manual.RDS")
  
  ## N = Number of series
  N <- 4
  
  # Tidy Up
  rm(mydata)
  gc(TRUE)

}
##====================== Data Load ============================##



##====================== Plot the Data =========================##
if (F) {
  ymin <- -15
  ymax <- 15
  par(mar = c(2.5, 2.5, 2, 1))
  
  ptitle <- "ANZ"
  plot(dates,e_anz,"l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax)) # Dates on x-axis
  ptitle <- "CBA"
  plot(dates,e_cba,"l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  ptitle <- "NAB"
  plot(dates,e_nab,"l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
  ptitle <- "WBC"
  plot(dates,e_wbc,"l",main = ptitle,ylab="",xlab="",font.main=1, ylim=c(ymin, ymax))
}
################# ABS(ret) and g PLOTS ########################
if(F){
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
if(F){
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
if(F){
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
#### Step 1:  Begin the multivariate estimation               ####
##==============================================================##

#Build the matrix of data:
e <- cbind(e_anz,e_cba,e_nab,e_wbc)

z_anz <- e_anz/sqrt(MTV_anz$Estimated$tv@g * MTV_anz$Estimated$garch@h)
z_cba <- e_cba/sqrt(MTV_cba$Estimated$tv@g * MTV_cba$Estimated$garch@h)
z_nab <- e_nab/sqrt(MTV_nab$Estimated$tv@g * MTV_nab$Estimated$garch@h)
z_wbc <- e_wbc/sqrt(MTV_wbc$Estimated$tv@g * MTV_wbc$Estimated$garch@h)

z <- cbind(z_anz,z_cba,z_nab,z_wbc)
Tobs <- NROW(z)

STCC <- stcc(Tobs)
STCC$P1 <- matrix(0.4,N,N)  # Const Corr = 0.1
diag(STCC$P1) <- 1
STCC$P2 <- matrix(0.7,N,N)  # Const Corr = 0.9
diag(STCC$P2) <- 1
STCC$TRpars <- c(2,0.5)     # (speed,location)
STCC$shape <- TVshape$single

## Modify start pars to improve estimates
STCC$P1[4,1] <- 0.55
STCC$P1[1,4] <- 0.55

nr.pars <- 14
STCC$optimcontrol$ndeps <- rep(1e-5,nr.pars)

timestamp()
STCC <- EstimateSTCC(z,STCC,calcHess = TRUE)
timestamp()

seP1 <- round(matrix(myUnVecl(STCC$Estimated$stderr[1:(N*(N-1)/2)]),N,N),4)
seP2 <- round(matrix(myUnVecl(STCC$Estimated$stderr[(N*(N-1)/2+1):(N*(N-1))]),N,N),4)
seTR <- round(tail(STCC$Estimated$stderr,2),3)

round(STCC$Estimated$P1,4)
seP1

round(STCC$Estimated$P2,4)
seP2

STCC$Estimated$TRpars
seTR

STCC$Estimated$value

saveRDS(STCC,"Output/STCC_Estimated.RDS")

## The end  ####

Pt <- STCC$Estimated$Pt

plot(Pt[,1],type = 'l',ylim=c(0.2,0.9),panel.first = grid(),main="Window:250,data:z")
lines(Pt[,2],type = 'l',col="red")
lines(Pt[,3],type = 'l',col="green")
lines(Pt[,4],type = 'l',col="blue")
lines(Pt[,5],type = 'l',col="darkgrey")
lines(Pt[,6],type = 'l',col="yellow",lwd=2)



## Calc local correlation using rolling window ####

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

Pt_rollingWin <- calculate_local_corr(z,250)
Pt <- Pt_rollingWin
# Re-run the charts using the code above




##==============================================================##
##                            THE END
##==============================================================##

