rm(list=ls()) 
##########  conditional independece test   ###############
library(MASS)
library(parallel)
library(glmnet)
library(cdcsis)
library(GeneralisedCovarianceMeasure)
library(devtools)
library(remotes)
#remotes::install_github("ericstrobl/RCIT", force = TRUE)
library(RCIT)


#source('D:/0-xnc-du-dr/CIN/all-functions-CDI.R')
source('all-functions-CDI-p.R')

#Rcpp::sourceCpp('D:/0-xnc-du-dr/0-my_codes/cv_hat_max_arma_r.cpp')
Rcpp::sourceCpp('cv_hat_max_arma_r.cpp')
###### data

#dat2 <- read.csv("D:/0-xnc-du-dr/data/GICS/X_Mat.csv")

alpha <- 0.05





Gauss_CI <- function(xx, yy, zz, alpha){
  x <- xx
  y <- yy
  z <- zz
  dn <- dim(x)[1]
  #if(dim(dat$x)[1] <= dn){
  #  dn= dim(dat$x)[1]
  #  x <- dat$x[1:dn, ]
  #  y <- dat$y[1:dn, ]
  #  z <- dat$z[1:dn, ]
  #}else{
  #  l <- dim(dat$x)[1]
  #  ind_c <- sample(l, dn, replace = F, prob=rep(1/l, l))
  #  x <- dat$x[ind_c, ]
  #  y <- dat$y[ind_c, ]
  #  z <- dat$z[ind_c, ]
  #}
  n <- dim(x)[1]
  p <- dim(x)[2]
  q <- dim(y)[2]
  m <- dim(z)[2]
  
  time_m <- rep(0, 8)
  Un_m <- rep(0, 8)
  ##################
  t1 <- Sys.time()
  z_norm <- data_transform(z)$x_trs_my
  x_norm <- data_transform(x)$x_trs_my
  y_norm <- data_transform(y)$x_trs_my
  if(dim(x)[2]==1){
    Epsilon <- matrix(lm(x_norm~z_norm-1)$residuals, n)
    Delta <- matrix(lm(y_norm~z_norm-1)$residuals, n)
  }else{
    Epsilon <- apply(x_norm, 2, fun_penalty, z_norm)
    Delta   <- apply(y_norm, 2, fun_penalty, z_norm)
  }
  #########
  #Epsi_norm <- data_transform(Epsilon)$x_trs_my
  #Delt_norm <- data_transform(Delta)$x_trs_my
  Epsi_norm <- Epsilon
  Delt_norm <- Delta 
  N <- 5000*n
  M <- NULL
  #M <- matrix(rnorm(N, 0, 1), 5000)
  M <- matrix(sample(c(-1,1), N,  replace = T, prob = rep(0.5, 2)), 5000)
  #############
  res <- cv_arma(Epsi_norm, Delt_norm, M, alpha)
  #scv_hat <- res$cv_est
  ##############test statistics
  sH_t <- sort(abs(1/n * apply(res$gamma, 2, sum)),decreasing = TRUE)
  #sH_t <- sort(abs(1/n * t(x_norm) %*% (y_norm)), decreasing = TRUE)
  sH_test <- rep(0, 3)
  sH_test[1] <- n^{1/2} * sH_t[1]
  sH_test[2] <- n^{1/2} * sum(sH_t[1 : 3])
  sH_test[3] <- n^{1/2} * sum(sH_t[1 : 5])
  tes_s1 <- res$tess[, 1]
  tes_s2 <- apply(res$tess[, 1:3], 1, sum)
  tes_s3 <- apply(res$tess[, 1:5], 1, sum)
  gauspVal <- rep(0, 3)
  gauspVal[1] <- mean(tes_s1 > sH_test[1]) 
  gauspVal[2] <- mean(tes_s2 > sH_test[2]) 
  gauspVal[3] <- mean(tes_s3 > sH_test[3])
  t2 <- Sys.time()
  time_m[6:8] <- t2-t1
  Un_m[6:8] <- attr(t2-t1,"unit")
  
  ###########zhu PCD
  t1 <- Sys.time()
  pdcv_inv <- try(pro_Val(x, y, z, alpha)$pval,  silent=TRUE)
  t2 <- Sys.time()
  if ('try-error' %in% class(pdcv_inv)) {
    PDCv <- 8888 
    time_m[1] <- 8888
  }else{ 
    PDCv <- pdcv_inv
    time_m[1] <- t2-t1
  }
  Un_m[1] <- attr(t2-t1,"unit")
  
  ######RCIT
  t1 <- Sys.time()
  rcitp <- RCIT::RCIT(x, y, z)$p
  ritv <-  rcitp
  t2 <- Sys.time()
  time_m[2] <- t2-t1
  Un_m[2] <- attr(t2-t1,"unit")
  ######RcoT
  t1 <- Sys.time()
  rcotp <- RCIT::RCoT(x, y, z)$p
  Rcotv <- rcotp
  t2 <- Sys.time()
  time_m[3] <- t2-t1
  Un_m[3] <- attr(t2-t1,"unit")
  
  ######
  ### cdsis
  t1 <- Sys.time()
  cdorPval <- cdcsis::cdcov.test(x, y, z)$p.value
  cdor <-  cdorPval
  t2 <- Sys.time()
  time_m[4] <- t2-t1
  Un_m[4] <- attr(t2-t1,"unit")
  ###1: reject    generalised
  t1 <- Sys.time()
  gcm <- as.numeric(GeneralisedCovarianceMeasure::gcm.test(x, y, z , alpha = 0.05, regr.method="kernel.ridge")$p.value)
  t2 <- Sys.time()
  time_m[5] <- t2-t1
  Un_m[5] <- attr(t2-t1,"unit")
  resn <- c(PDCv, ritv, Rcotv, cdor, gcm, gauspVal)
  result <- list(resn = resn, time = time_m, Un_m = Un_m, dn=dn, p=p, q=q, m=m)
  return(result)
}


name <- c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
          'Ind', 'IT', 'Mat', 'RE', 'Uti')

#name <- c('CS', 'CD', 'CSt')
com_my <- combn(name, 2)


start <- Sys.time()


rep.num <- dim(com_my)[2]
mycores <- 20
res_gauss <- mclapply(1:rep.num, function(rep.sim){
  set.seed(12345)
  i <- rep.sim
  res  <- rep(0, 13)
  coln <- paste0(com_my[1,i], "-", com_my[2,i])
  a1 <- paste0("02return_", com_my[1,i], ".csv")
  dat1 <- read.csv(a1)
  a2 <- paste0("02return_", com_my[2,i], ".csv")
  dat2 <- read.csv(a2)
  xx <- NULL;  yy <- NULL
  xx <- as.matrix(dat1)[,-1]
  yy <- as.matrix(dat2)[,-1]
  ######zz
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
  res_result <- Gauss_CI(xx, yy, zz, alpha)
  res[1] <- coln
  res[2] <- res_result$dn
  dim_m <- matrix(c(res_result$p, res_result$q, res_result$m), 3, 1)
  res[3:13] <- rbind(dim_m, matrix(res_result$resn,8))
  print(rep.sim)
  #write.csv(res_result, paste0("para_res", "_rep.sim", rep.sim, ".csv"))
  return(res)
}
, mc.cores = mycores)

res_1 <- matrix(0, 12, rep.num)
coln <- rep(0, rep.num)
for (i in 1:rep.num) {
  coln[i] <- res_gauss[[i]][1]
  res_1[,i] <- as.numeric(res_gauss[[i]][2:13])
}
rownames(res_1) <- c("ndim", "xdim", "ydim", "zdim" ,
                   "PCD", "PCIT", "PCOT", "CDSIS", "GCD", 
                   "gauss-1",   "gauss-3",   "gauss-5")

colnames(res_1) <- coln

print(t(res_1))


write.csv(t(res_1), paste0("p-cin-other-rescpval",".csv"))





