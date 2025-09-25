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


source('cind-functions.R') 
source('PCIND.R')

Rcpp::sourceCpp('cv_hat_max_arma.cpp')

alpha <- 0.05

Gauss_CI_data <- function(xx, yy, zz, alpha){
  x <- xx
  y <- yy
  z <- zz
  dn <- dim(x)[1]
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  q <- dim(y)[2]
  m <- dim(z)[2]
  sA <- c(0,0,0)
  PDCv <- ritv <- Rcotv <- cdor <- gcm <- 0
  
  ############# the proposed conditional independence test based on linear regressions ############
  seed_my = 12345
  res_my <- Cind_Gtest(x, y, z, alpha, seed_my, "all")
  sA[1] <- res_my$Gaussian$p_value      # Gaussian multiplier
  sA[2] <- res_my$Mammen$p_value        # Mammen's multiplier 
  sA[3] <- res_my$Rademacher$p_value    # Radamacher multiplier
  ##################
  
  ####################### PCD ########################
  set.seed(12345)
  pdcv_inv <- try(pro_Val(x, y, z, alpha)$pval,  silent=TRUE)
  t2 <- Sys.time()
  if ('try-error' %in% class(pdcv_inv)) {
    PDCv <- 8888        # indicates failure; downstream code must check and skip
  }else{
    PDCv <- pdcv_inv
  }
  
  ##################### RCIT ########################
  set.seed(12345)
  rcitp <- RCIT::RCIT(x, y, z,seed = 12345)$p
  ritv <-  rcitp
  
  ##################### RCoT ############
  set.seed(12345)
  rcotp <- RCIT::RCoT(x, y, z,seed = 12345)$p
  Rcotv <- rcotp
  
  
  ############### cdCov ##################### 
  set.seed(12345)
  cdorPval <- cdcsis::cdcov.test(x, y, z)$p.value
  cdor <-  cdorPval
  
  ##################### GCM #####################
  set.seed(12345)
  gcm <- as.numeric(GeneralisedCovarianceMeasure::gcm.test(x, y, z , alpha = 0.05, regr.method="kernel.ridge")$p.value)
  
  resn <- c(sA, gcm, PDCv, ritv, Rcotv, cdor)
  return(resn)
}


#-------------------- before COVID-19 period: 1 January 2016 - 31 December 2018 ----------------------------

name <- c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
          'Ind', 'IT', 'Mat', 'RE', 'Uti')

 
com_my <- combn(name, 2)


start <- Sys.time()


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
  ###### read z
  ind1 <- which(name==com_my[1,i])
  ind2 <- which(name==com_my[2,i])
  name_z <- name[-c(ind1, ind2)]
  a3 <- paste0("68return_", name_z[1], ".csv")
  dat_z <- as.matrix(read.csv(a3))[, -1]
  if(length(name_z) > 1){
    for (j in 2:length(name_z)) {
      az <- paste0("68return_", name_z[j], ".csv")
      dat_e <- read.csv(az)
      dat_z <- cbind(dat_z, as.matrix(dat_e)[,-1])
    }
  }
  zz <- dat_z
  # implement the proposed method (CI-Lasso) and competing methods to test the conditional independence between x and y given z
  res_result <- Gauss_CI_data(xx, yy, zz, alpha)
  res <-  matrix(res_result,8)
  print(rep.sim)
  return(res)
}
, mc.cores = mycores)

res_1 <- matrix(0, 8, rep.num)
for (i in 1:rep.num) {
  res_1[,i] <- as.numeric(res_gauss[[i]])
}


Name1 <- as.character(com_my[1, ])
Name2 <- as.character(com_my[2, ])

res_new <- t(res_1)  

colnames(res_new) <- c("Gaussian","Mammen","Rademacher","GCM", "PCD","PCIT","PCoT","cdCov")

res_all <- data.frame(
  name1 = Name1,
  name2 = Name2,
  res_new,
  check.names = FALSE,            
  stringsAsFactors = FALSE        
)


write.csv(res_all, paste0("68pval",".csv"))

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
  ###### read z
  ind1 <- which(name==com_my[1,i])
  ind2 <- which(name==com_my[2,i])
  name_z <- name[-c(ind1, ind2)]
  a3 <- paste0("02return_", name_z[1], ".csv")
  dat_z <- as.matrix(read.csv(a3))[, -1]
  if(length(name_z) > 1){
    for (j in 2:length(name_z)) {
      az <- paste0("02return_", name_z[j], ".csv")
      dat_e <- read.csv(az)
      dat_z <- cbind(dat_z, as.matrix(dat_e)[,-1])
    }
  }
  zz <- dat_z
  # implement the proposed method (CI-Lasso) and competing methods to test the conditional independence between x and y given z
  res_result <- Gauss_CI_data(xx, yy, zz, alpha)
  res <-  matrix(res_result,8)
  print(rep.sim)
  return(res)
}
, mc.cores = mycores)

res_1 <- matrix(0, 8, rep.num)
for (i in 1:rep.num) {
  res_1[,i] <- as.numeric(res_gauss[[i]])
}


Name1 <- as.character(com_my[1, ])
Name2 <- as.character(com_my[2, ])

res_new <- t(res_1)  

colnames(res_new) <- c("Gaussian","Mammen","Rademacher","GCM", "PCD","PCIT","PCoT","cdCov")

res_all <- data.frame(
  name1 = Name1,
  name2 = Name2,
  res_new,
  check.names = FALSE,            
  stringsAsFactors = FALSE        
)


write.csv(res_all, paste0("02pval",".csv"))
