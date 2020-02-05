
# setwd("C:/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project10_Eigentest")        #Anna's work PC - GOOGLE DRIVE
#setwd("E:/WORK/GoogleDrive/MyDocs/Annastiina/Topics/MTVGJR_MGARCH/Project10_Eigentest")   #Anna's Laptop - GOOGLE DRIVE
setwd("D:/OneDrive/Documents/Annastiina/Topics/MTVGJR_MGARCH/Project2_Building")    # Glen's home-office PC

# Results from mid-March, 2018
# Used this file to generate: AUS_4Banks_MultiVar_Estimation_v0.1.r
FourBanksResults1 <- readRDS("Results/MTVGJR_Estimated.RDS")
FourBanksResults2 <- readRDS("Results/MTVGJR_Estimated_vartarget_OFF.RDS")
FourBanksResults3 <- readRDS("Results/MTVGJR_Estimated_vartarget_OFF_2018.RDS")


# Results from late-March, 2018
# Used this file to generate: AUS_4Banks_MultiVar_Estimation_v1.0.r
FourBanksResults4 <- readRDS("Results/TVGJR_MultivarEstimation_1.RDS")
FourBanksResults5 <- readRDS("Results/TVGJR_MultivarEstimation_2.RDS")
FourBanksResults6 <- readRDS("Results/TVGJR_MultivarEstimation_2018_1.RDS")
FourBanksResults7 <- readRDS("Results/TVGJR_MultivarEstimation_2018_2.RDS")
FourBanksResults8 <- readRDS("Results/TVGJR_MultivarEstimation_2018_3.RDS")

# Analysis for paper - plots etc. + attempt at Std Error Calcs: AUS_4Banks_MultiVar_Estimation_v2.0.r

# More Analysis + attempt at Std Error Calcs: AUS_4Banks_MultiVar_Estimation_2018.r


View(FourBanksResults1)