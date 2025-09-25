rm(list=ls()) 
##########  conditional independece test   ###############
library(MASS)
library(parallel)
library(glmnet)
library(cdcsis)
library(GeneralisedCovarianceMeasure)
library(remotes)
#remotes::install_github("ericstrobl/RCIT", force = TRUE)
library(RCIT)


source('PIND.R')

Rcpp::sourceCpp('cv_hat_max_arma.cpp')

alpha <- 0.05




#-------------------- before COVID-19 period: 1 January 2016 - 31 December 2018 ----------------------------

name <- c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
          'Ind', 'IT', 'Mat', 'RE', 'Uti')


com_my <- combn(name, 2)

rep.num <- dim(com_my)[2]
mycores <- 55
res_gauss <- mclapply(1:rep.num, function(rep.sim){
  set.seed(12345)
  # read x and y
  i <- rep.sim
  res  <- rep(0, 8)
  coln <- paste0(com_my[1,i], "-", com_my[2,i])
  a1 <- paste0("68return_", com_my[1,i], ".csv")
  dat1 <- read.csv(a1)
  a2 <- paste0("68return_", com_my[2,i], ".csv")
  dat2 <- read.csv(a2)
  xx <- NULL;  yy <- NULL
  xx <- as.matrix(dat1)[,-1]
  yy <- as.matrix(dat2)[,-1]
  
  # implement the proposed method test the independence between x and y

  sA <- c(0,0,0)
  seed_my = 12345
  res_my <- Ind_Gtest(xx, yy, alpha, seed_my, "all")
  sA[1] <- res_my$Gaussian$p_value      # Gaussian multiplier
  sA[2] <- res_my$Mammen$p_value        # Mammen's multiplier 
  sA[3] <- res_my$Rademacher$p_value    # Radamacher multiplier
  print(rep.sim)
  return(sA)
}
, mc.cores = mycores)

res_1 <- matrix(0, 3, rep.num)
for (i in 1:rep.num) {
  res_1[,i] <- as.numeric(res_gauss[[i]])
}


Name1 <- as.character(com_my[1, ])
Name2 <- as.character(com_my[2, ])

res_new <- t(res_1)  

colnames(res_new) <- c("Gaussian","Mammen","Rademacher")

res_all <- data.frame(
  name1 = Name1,
  name2 = Name2,
  res_new,
  check.names = FALSE,            
  stringsAsFactors = FALSE        
)


write.csv(res_all, paste0("i68pval",".csv"))

#------------------------- during/after COVID-19 period: 1 January 2020 - 31 December 2022 ---------------------------------------------

rep.num <- dim(com_my)[2]
mycores <- 55
res_gauss <- mclapply(1:rep.num, function(rep.sim){
  set.seed(12345)
  # read x and y
  i <- rep.sim
  res  <- rep(0, 8)
  coln <- paste0(com_my[1,i], "-", com_my[2,i])
  a1 <- paste0("02return_", com_my[1,i], ".csv")
  dat1 <- read.csv(a1)
  a2 <- paste0("02return_", com_my[2,i], ".csv")
  dat2 <- read.csv(a2)
  xx <- NULL;  yy <- NULL
  xx <- as.matrix(dat1)[,-1]
  yy <- as.matrix(dat2)[,-1]
  sA <- c(0,0,0)
  
  # implement the proposed method test the independence between x and y
  
  seed_my = 12345
  res_my <- Ind_Gtest(xx, yy, alpha, seed_my, "all")
  sA[1] <- res_my$Gaussian$p_value      # Gaussian multiplier
  sA[2] <- res_my$Mammen$p_value        # Mammen's multiplier 
  sA[3] <- res_my$Rademacher$p_value    # Radamacher multiplier
  print(rep.sim)
  return(sA)
}
, mc.cores = mycores)

res_1 <- matrix(0, 3, rep.num)
for (i in 1:rep.num) {
  res_1[,i] <- as.numeric(res_gauss[[i]])
}


Name1 <- as.character(com_my[1, ])
Name2 <- as.character(com_my[2, ])

res_new <- t(res_1)  

colnames(res_new) <- c("Gaussian","Mammen","Rademacher")

res_all <- data.frame(
  name1 = Name1,
  name2 = Name2,
  res_new,
  check.names = FALSE,            
  stringsAsFactors = FALSE        
)


write.csv(res_all, paste0("i02pval",".csv"))
