rm(list=ls()) 
##########  independece test   ###############
#setwd('D:/0-xnc-du-dr/code')
library(MASS)
library(HHG)
library(energy)
library(dHSIC)
library(clue)
library(adagio)
library(pracma)
library(randtoolbox)
library(EDMeasure)
library(mvtnorm)
library(purrr)
library(parallel)
require(randtoolbox, quietly = T)


#source('D:/0-xnc-du-dr/1-paper/IN-result/all-functions.R')
source('all-functions.R')

#Rcpp::sourceCpp('D:/0-xnc-du-dr/1-paper/IN-result/cv_hat_max_arma_r.cpp')
Rcpp::sourceCpp('cv_hat_max_arma_r.cpp')
###### data

#dat2 <- read.csv("D:/0-xnc-du-dr/data/GICS/X_Mat.csv")


alpha <- 0.05

Gauss_T <- function(xx, yy, zz, alpha){
  x <- xx
  y <- yy
  dn <- dim(x)[1]
  n <- dim(x)[1]
  p <- dim(x)[2]
  q <- dim(y)[2]
  sA <- rep(0, 3)
  sB <- rep(0, 3)
  Jd <- rep(0, 3)
  time_m <- rep(0, 14)
  Un_m <- NULL
  BstrpTimes <- 500
  Bstrp <- matrix(0,BstrpTimes,1)
  Bs_han <- matrix(0,BstrpTimes,1)
  Bs_Mr <- matrix(0,BstrpTimes,1)
  t1 <- Sys.time()
  x_norm <- data_transform(x)$x_trs_my
  y_norm <- data_transform(y)$x_trs_my
  N <- 5000*n
  M <- NULL
  M <- matrix(sample(c(-1,1), N,  replace = T, prob = rep(0.5, 2)), 5000)
  #############
  res <- cv_arma(x_norm, y_norm, M, alpha)
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
  #gauspVal[2] <- mean(tes_s2 > sH_test[2])
  gauspVal[2] <- sH_test[1]
  #gauspVal[3] <- mean(tes_s3 > sH_test[3])
  t2 <- Sys.time()
  time_m[12:14] <- t2-t1
  Un_m[12:14] <- attr(t2-t1,"unit")
  
  
  ########
  t1 <- Sys.time()
  #Indep_Index =My_two_MULTI_U(x,y);
  #for (jj in 1:BstrpTimes){
  #  O=order(runif(n,0,1));
  #  YBstrp = y[O,];
  #  Bstrp[jj] = My_two_MULTI_U(x,YBstrp);
  #}
  pVal = NA #mean(Bstrp>=Indep_Index);
  t2 <- Sys.time()
  time_m[1] <- t2-t1
  Un_m[1] <- attr(t2-t1,"unit")
  
  #######
  #t1 <- Sys.time()
  #Dx = as.matrix(dist((x)), diag = TRUE, upper = TRUE)
  #Dy = as.matrix(dist((y)), diag = TRUE, upper = TRUE)
  #hhg = hhg.test(Dx, Dy, nr.perm =BstrpTimes)
  drpVal = NA #hhg$perm.pval.hhg.sc;
  t2 <- Sys.time()
  time_m[2] <- t2-t1
  Un_m[2] <- attr(t2-t1,"unit")
  
  ######
  t1 <- Sys.time()
  dcortpVal =dcorT.test(x,y)$p.value;
  t2 <- Sys.time()
  time_m[3] <- t2-t1
  Un_m[3] <- attr(t2-t1,"unit")
  
  ###########dhsic
  t1 <- Sys.time()
  dhsicpVal <- dhsic.test(x, y)$p.value
  t2 <- Sys.time()
  time_m[4] <- t2-t1
  Un_m[4] <- attr(t2-t1,"unit")
  
  
  z <- cbind(x,y)
  d <- dim(z)[2]
  
  
  #####################jd
  Mp=500
  sta13 <- function(x, Mp, n, d, alpha){
    res <- rep(0,3)
    Ustat=TRUE
    cc=1
    if(Ustat) F <- Jdcov.sq.U else F <- Jdcov.sq.V
    X <- split(x, rep(1:ncol(x), each = nrow(x)))
    x.r <- x
    for(j in 1:d) { f.cdf <- ecdf(x[,j]); x.r[,j] <- f.cdf(x[,j]) }
    
    stat <- F(x,cc)
    stat.r <- F(x.r,cc)
    stat.s <- Jdcov.sq.US(x,cc)
    stat.p <- stat.pr <- stat.ps <- rep(0, Mp)
    for(i in 1:Mp)
    {
      x.pr <- x.p <- x
      for(j in 1:d) x.p[,j] <- x[sample(1:n,n,replace=TRUE),j]
      X.p <- split(x.p, rep(1:ncol(x.p), each = nrow(x.p)))
      for(j in 1:d) { f.cdf <- ecdf(x.p[,j]); x.pr[,j] <- f.cdf(x.p[,j]) }
      stat.p[i] <- F(x.p,cc)
      stat.pr[i] <- F(x.pr,cc)
      stat.ps[i] <- Jdcov.sq.US(x.p,cc)
      
    }
    
    #crit   <- quantile(stat.p, c(1-alpha))
    #crit.r <- quantile(stat.pr, c(1-alpha))
    #crit.s <- quantile(stat.ps, c(1-alpha))
    #res[1]  <- ifelse(stat > crit, 1, 0)
    #res[2]  <- ifelse(stat.r > crit.r, 1, 0)
    #res[3]  <- ifelse(stat.s > crit.s, 1, 0)
    res[1]  <- mean(stat.p > stat)
    res[2]  <- mean(stat.pr > stat.r)
    res[3]  <- mean(stat.ps > stat.s)
    return(res)
  }
  t1 <- Sys.time()
  time_m[4] <- t2-t1
  #jdcov_inv <- try(sta13(z, Mp, n, d, alpha),  silent=TRUE) 
  t2 <- Sys.time()
  #if ('try-error' %in% class(jdcov_inv)) {
    Jd <- rep(NA,3)
    time_m[5:7] <- NA
    #} else{ 
    #  Jd <- jdcov_inv
    #  time_m[5:7] <- t2-t1
    # }
  Un_m[5:7] <- attr(t2-t1,"unit")
  ########TMT
  sta4 <- function(x, Mp, n, d, alpha){
    stat4 <- 0
    stat.p4 <- rep(0, Mp)
    
    for(j in 1:(d-1)) stat4 <- stat4 + dcovU(x[,j],x[,(j+1):d])^2
    ##
    for(i in 1:Mp)
    {
      x.pr <- x.p <- x
      for(j in 1:d) x.p[,j] <- x[sample(1:n,n,replace=TRUE),j]
      X.p <- split(x.p, rep(1:ncol(x.p), each = nrow(x.p)))
      
      for(j in 1:(d-1)) stat.p4[i] <- stat.p4[i] + dcovU(x.p[,j],x.p[,(j+1):d])^2
    }
    #crit4 <- quantile(stat.p4, c(1-alpha))
    #TMT <- ifelse(stat4 > crit4, 1, 0)
    TMT <- mean(stat.p4 > stat4)
    return(TMT)
  }
  t1 <- Sys.time()
  #tmt.inv <- try(sta4(z, Mp, n, d, alpha), silent=TRUE)
  t2 <- Sys.time()
  #if ('try-error' %in% class(tmt.inv)) {
    TMT <- NA 
    time_m[8] <- NA
    #} else{ 
    #TMT <- tmt.inv
    #time_m[8] <- t2-t1
    #}
  Un_m[8] <- attr(t2-t1,"unit")
  
  ####mdm
  t1 <- Sys.time()
  #m.inv <- try(mdm_test(z, c(p,q))$pval, silent=TRUE)  
  t2 <- Sys.time()
  #if ('try-error' %in% class(m.inv)) {
    dmdmpVal <- NA 
    time_m[9] <- NA
    #} else{ 
    #  dmdmpVal <- m.inv 
    #  time_m[9] <- t2-t1
    #}
  Un_m[9] <- attr(t2-t1,"unit")
  
  ######### Han Hallin
  t1 <- Sys.time()
  #had.inv <- try({ 
  #   H_dcovt <- discov.Hallin.stat(x,y)
    #   for (j in 1:BstrpTimes){
      #    O1=order(runif(n,0,1))
      #    YBstrp1 = y[O1,]
      #    Bs_han[j]  = discov.Hallin.stat(x,YBstrp1)
      #  } 
    #} , silent=TRUE)
  t2 <- Sys.time()
  #if ('try-error' %in% class(had.inv)) {
    HaVal <- NA 
    time_m[10] <- NA
    #}else{ 
    # HaVal <- mean(Bs_han >= H_dcovt)
    # time_m[10] <- t2-t1
    #}
  Un_m[10] <- attr(t2-t1,"unit")
  
  ######Mrcov
  t1 <- Sys.time()
  #mrn.inv <- try( { 
  #   Mrcov <- computestatisticrdcov(x, y)
    #   for (j in 1:BstrpTimes){
      #     O2=order(runif(n,0,1));
      #     YBstrp2 = y[O2,];
      #     Bs_Mr[j]  = computestatisticrdcov(x, YBstrp2);
      #   }
    #} , silent=TRUE)
#t2 <- Sys.time()
  #if ('try-error' %in% class(mrn.inv)) {
    MrVal <- NA 
    time_m[11] <- NA
    # }else{ 
    #  MrVal <- mean(Bs_Mr >= Mrcov)
    #  time_m[11] <- t2-t1
    #}
  Un_m[11] <- attr(t2-t1,"unit")
  
  resn <- c(pVal, drpVal, dcortpVal, dhsicpVal,
            Jd, TMT,  dmdmpVal, HaVal, MrVal, gauspVal)
  result <- list(resn = resn, time = time_m, Un_m = Un_m, dn=dn, p=p, q=q)
  return(result)
}


name <- c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
          'Ind', 'IT', 'Mat', 'RE', 'Uti')

#name <- c('CS', 'CD')
com_my <- combn(name, 2)


start <- Sys.time()
rep.num <- dim(com_my)[2]
mycores <- 1
res_gauss <- mclapply(1:rep.num, function(rep.sim){
set.seed(12345)
  i <- rep.sim
  res  <- rep(0, 18)
  coln <- paste0(com_my[1,i], "-", com_my[2,i])
  a1 <- paste0("02return_", com_my[1,i], ".csv")
  dat1 <- read.csv(a1)
  a2 <- paste0("02return_", com_my[2,i], ".csv")
  dat2 <- read.csv(a2)
  xx <- NULL;  yy <- NULL
  xx <- as.matrix(dat1)[,-1]
  yy <- as.matrix(dat2)[,-1]
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
  res_result <- Gauss_T(xx, yy, zz, alpha)
  res[1] <- coln
  res[2] <- res_result$dn
  dim_m <- matrix(c(res_result$p, res_result$q), 2, 1)
  res[3:18] <- rbind(dim_m, matrix(res_result$resn,14))
  print(rep.sim)
  #write.csv(res_result, paste0("para_res", "_rep.sim", rep.sim, ".csv"))
  return(res)
  }
, mc.cores = mycores)

res_1 <- matrix(0, 17, rep.num)
coln <- rep(0, rep.num)
for (i in 1:rep.num) {
  coln[i] <- res_gauss[[i]][1]
  res_1[,i] <- as.numeric(res_gauss[[i]][2:18])
}
rownames(res_1) <- c("ndim", "xdim", "ydim",  
                   "pcov", "drcov", "dcor", "hsic", 
                   "Jdcov", "Jdcov_r", "Jdcov_s", 
                   "TMT", "mdm", "Hallin", "Mrank",
                   "gauss-1",   "gauss-3",   "gauss-5")

colnames(res_1) <- coln

print(t(res_1))


write.csv(t(res_1), paste0("1respv",".csv"))
#write.csv(t(res_1), "D:/0-xnc-du-dr/1-paper/sta3.csv")  

