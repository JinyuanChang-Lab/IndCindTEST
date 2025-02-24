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


#source('D:/0-xnc-du-dr/0-my_codes/all-functions.R')
source('all-functions.R')

#Rcpp::sourceCpp('D:/0-xnc-du-dr/1code/0307/cpp/cv_hat_max_arma.cpp')
Rcpp::sourceCpp('cv_hat_max_arma.cpp')
###### data

rmammen <- function(n,
                    construct = c("normal-2", "normal-1", "two-point mass")){
  construct <- match.arg(construct)
  
  if (length(n) > 1) n <- length(n)
  
  if (construct == "two-point mass"){
    .vals <- c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2)
    .probs <- rev(abs(.vals)/sqrt(5))
    
    sample(.vals, size = n, replace = TRUE, prob = .probs)
  } else if (construct == "normal-1"){
    .v <- rnorm(n)
    
    (.v/sqrt(2)) + (0.5*((.v*.v) - 1))
  } else {
    .delta1 <- sqrt((3/4) + (sqrt(17)/12))
    .delta2 <- sqrt((3/4) - (sqrt(17)/12))
    .v <- matrix(rnorm(2*n), nrow=n, ncol=2)
    
    (.delta1 + (.v[,1]/sqrt(2)))*(.delta2 + (.v[,2]/sqrt(2))) -
      (.delta1*.delta2)
  }
}



Data1 <- function(n, p, q, dep){
  if(dep==0){
    x <- 0.2 * matrix(rcauchy(n*p, 0, 1), n, p) 
    y <- 0.2 * matrix(rcauchy(n*p, 0, 1), n, p)
  }else{
    if(dep==1){d <- p/20
    }else if(dep==2){d <- p/10
    }else{d <- p/5}
    x <- 0.2 * matrix(rcauchy(n*p, 0, 1), n, p)
    y <- 0.2 * matrix(rcauchy(n*p, 0, 1), n, p)
    prep <- rnorm(n,0,1)
    A <- matrix(rep(prep, d), n)
    x[,1:d] <- x[,1:d] + A
    y[,1:d] <- y[,1:d] + A
  }
  dat <- list(x=x, y=y)
  return(dat)  
}
alpha <- 0.05
n <- 50
p <- 100
q <- 100
###1,p/20; 2,p/10; 3,p/5
dep <- 0

set.seed(12345)

Gauss_T <- function(n, p, q, dep, alpha){
  sA <- rep(0, 12)
  sB <- rep(0, 12)
  Jd <- rep(0, 3)
  time_m <- rep(0, 35)
  Un_m <- rep(0, 35)
  BstrpTimes <- 500
  Bstrp <- matrix(0,BstrpTimes,1)
  Bs_han <- matrix(0,BstrpTimes,1)
  Bs_Mr <- matrix(0,BstrpTimes,1)
  dat <- Data1(n, p, q, dep)
  x <- dat$x
  y <- dat$y
  t1 <- Sys.time()
  x_norm <- data_transform(x)$x_trs_my
  y_norm <- data_transform(y)$x_trs_my
  N <- 5000*n
  M <- NULL
  M <- matrix(sample(c(-1,1), N,  replace = T, prob = rep(0.5, 2)), 5000)
  #############
  res <- cv_arma(x_norm, y_norm, M, alpha)$cv_est
  scv_hat <- res
  ##############test statistics
  sH_t <- sort(abs(1/n * t(x_norm) %*% (y_norm)), decreasing = TRUE)
  sH_test <- rep(0, 3)
  sH_test[1] <- n^{1/2} * sH_t[1]
  sA[1] <- ifelse(sH_test[1] > scv_hat[1], 1, 0)
  
  ###########Gaussian multipliers
  M_n <- matrix(rnorm(N, 0, 1), 5000)
  #############
  res_n <- cv_arma(x_norm, y_norm, M_n, alpha)$cv_est
  scv_hat_n <- res_n
  ##############test statistics
  sH_n <- sort(abs(1/n * t(x_norm) %*% (y_norm)), decreasing = TRUE)
  sH_test_n <- n^{1/2} * sH_n[1]
  sA[2] <- ifelse(sH_test_n > scv_hat_n[1], 1, 0)
  
  ###########rmammen  multipliers-"two-point mass"
  M_mt <- matrix(rmammen(N, "two-point mass"), 5000)
  #############
  res_mt <- cv_arma(x_norm, y_norm, M_mt, alpha)$cv_est
  scv_hat_mt <- res_mt
  ##############test statistics
  sH_mt <- sort(abs(1/n * t(x_norm) %*% (y_norm)), decreasing = TRUE)
  sH_test_mt <- n^{1/2} * sH_mt[1]
  sA[3] <- ifelse(sH_test_mt > scv_hat_mt[1], 1, 0)
  
  t2 <- Sys.time()
  time_m[12:35] <- t2-t1
  Un_m[12:35] <- attr(t2-t1,"unit")
  if(Un_m[12]=="mins" && p==400){
    time_m[12:35] <- time_m[12:35]*60
    Un_m[12:35] <- rep("secs", 24)
  }
  
  
  ###0.99
  scv_99 <- scv_hat * 0.99
  sA[4] <- ifelse(sH_test[1] > scv_99[1], 1, 0)
  sA[5] <- ifelse(sH_test[2] > scv_99[2], 1, 0)
  sA[6] <- ifelse(sH_test[3] > scv_99[3], 1, 0)
  ###0.98
  scv_98 <- scv_hat * 0.98
  sA[7] <- ifelse(sH_test[1] > scv_98[1], 1, 0)
  sA[8] <- ifelse(sH_test[2] > scv_98[2], 1, 0)
  sA[9] <- ifelse(sH_test[3] > scv_98[3], 1, 0)
  
  ###0.95
  scv_95 <- scv_hat * 0.95
  sA[10] <- ifelse(sH_test[1] > scv_95[1], 1, 0)
  sA[11] <- ifelse(sH_test[2] > scv_95[2], 1, 0)
  sA[12] <- ifelse(sH_test[3] > scv_95[3], 1, 0)
  
  
  ###########stand
  x_mean_n <- apply(x_norm, 2, mean)
  y_mean_n <- apply(y_norm, 2, mean)
  x_new_n <- x_norm - matrix(rep(x_mean_n, n), n, byrow=T) 
  y_new_n <- y_norm - matrix(rep(y_mean_n, n), n, byrow=T) 
  
  sH_ts <- sort(abs(1/n * t(x_new_n) %*% (y_new_n)), decreasing = TRUE)
  sH_tests <- rep(0, 3)
  sH_tests[1] <- n^{1/2} * sH_ts[1]
  sH_tests[2] <- n^{1/2} * sum(sH_ts[1 : 3])
  sH_tests[3] <- n^{1/2} * sum(sH_ts[1 : 5])
  ###########
  sB[1] <- ifelse(sH_tests[1] > scv_hat[1], 1, 0)
  sB[2] <- ifelse(sH_tests[2] > scv_hat[2], 1, 0)
  sB[3] <- ifelse(sH_tests[3] > scv_hat[3], 1, 0)
  ###0.99
  
  sB[4] <- ifelse(sH_tests[1] > scv_99[1], 1, 0)
  sB[5] <- ifelse(sH_tests[2] > scv_99[2], 1, 0)
  sB[6] <- ifelse(sH_tests[3] > scv_99[3], 1, 0)
  ###0.98
  
  sB[7] <- ifelse(sH_tests[1] > scv_98[1], 1, 0)
  sB[8] <- ifelse(sH_tests[2] > scv_98[2], 1, 0)
  sB[9] <- ifelse(sH_tests[3] > scv_98[3], 1, 0)
  
  ###0.95
  
  sB[10] <- ifelse(sH_tests[1] > scv_95[1], 1, 0)
  sB[11] <- ifelse(sH_tests[2] > scv_95[2], 1, 0)
  sB[12] <- ifelse(sH_tests[3] > scv_95[3], 1, 0)  
  
  ########
  t1 <- Sys.time()
  pcov_fun <- function(x,y,BstrpTimes){
    Bstrp <- matrix(0,BstrpTimes,1)
    Indep_Index =My_two_MULTI_U(x,y);
    for (jj in 1:BstrpTimes){
      O=order(runif(n,0,1));
      YBstrp = y[O,];
      Bstrp[jj]  =My_two_MULTI_U(x,YBstrp);
    }
    pVal = mean(Bstrp>=Indep_Index)
    return(pVal)
  }
  pcov_inv <- try(pcov_fun(x,y,BstrpTimes),  silent=TRUE) 
  t2 <- Sys.time()
  if ('try-error' %in% class(pcov_inv)) {
    pVal <- 8888 
    time_m[1] <- 8888
  }else{ 
    pVal = pcov_inv
    time_m[1] <- t2-t1
  }
  Un_m[1] <- attr(t2-t1,"unit")
  
  
  
  #######
  t1 <- Sys.time()
  dr_fun <- function(x,y, BstrpTimes){
    Dx = as.matrix(dist((x)), diag = TRUE, upper = TRUE)
    Dy = as.matrix(dist((y)), diag = TRUE, upper = TRUE)
    hhg = hhg.test(Dx, Dy, nr.perm =BstrpTimes)
    drpVal =hhg$perm.pval.hhg.sc;
    return(drpVal)
  } 
  dr_inv <- try(dr_fun(x,y, BstrpTimes),  silent=TRUE) 
  t2 <- Sys.time()
  if ('try-error' %in% class(dr_inv)) {
    drpVal <- 8888 
    time_m[2] <- 8888
  }else{ 
    drpVal <- ifelse(is.nan(dr_inv), 9999, dr_inv)
    time_m[2] <- t2-t1
  }
  Un_m[2] <- attr(t2-t1,"unit")
  
  ######
  t1 <- Sys.time()
  dcor_inv <- try(dcorT.test(x,y)$p.value,  silent=TRUE)
  t2 <- Sys.time()
  if ('try-error' %in% class(dcor_inv)) {
    dcortpVal <- 8888 
    time_m[3] <- 8888
  }else{ 
    dcortpVal <- ifelse(is.nan(dcor_inv), 9999, dcor_inv)
    time_m[3] <- t2-t1
  }
  Un_m[3] <- attr(t2-t1,"unit")
  
  ###########dhsic
  
  t1 <- Sys.time()
  dhsic_inv <- try(dhsic.test(x, y)$p.value,  silent=TRUE)
  t2 <- Sys.time()
  if ('try-error' %in% class(dhsic_inv)) {
    dhsicpVal <- 8888 
    time_m[4] <- 8888
  }else{ 
    dhsicpVal <- ifelse(is.nan(dhsic_inv)|is.na(dhsic_inv), 9999, dhsic_inv)
    time_m[4] <- t2-t1
  }
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
    
    #stat <- F(x,cc)
    stat.r <- F(x.r,cc)
    #stat.s <- Jdcov.sq.US(x,cc)
    #stat.p <- stat.pr <- stat.ps <- rep(0, Mp)
    for(i in 1:Mp)
    {
      x.pr <- x.p <- x
      for(j in 1:d) x.p[,j] <- x[sample(1:n,n,replace=TRUE),j]
      X.p <- split(x.p, rep(1:ncol(x.p), each = nrow(x.p)))
      for(j in 1:d) { f.cdf <- ecdf(x.p[,j]); x.pr[,j] <- f.cdf(x.p[,j]) }
      #stat.p[i] <- F(x.p,cc)
      stat.pr[i] <- F(x.pr,cc)
      #stat.ps[i] <- Jdcov.sq.US(x.p,cc)
      
    }
    
    #crit   <- quantile(stat.p, c(1-alpha))
    crit.r <- quantile(stat.pr, c(1-alpha))
    #crit.s <- quantile(stat.ps, c(1-alpha))
    #res[1]  <- ifelse(stat > crit, 1, 0)
    res[2]  <- ifelse(stat.r > crit.r, 1, 0)
    #res[3]  <- ifelse(stat.s > crit.s, 1, 0)
    return(res)
  }
  t1 <- Sys.time()
  jdcov_inv <- try(sta13(z, Mp, n, d, alpha),  silent=TRUE) 
  t2 <- Sys.time()
  if ('try-error' %in% class(jdcov_inv)) {
    Jd <- rep(8888, 3)
    time_m[5:7] <- 8888
  } else{ 
    Jd <- jdcov_inv
    time_m[5:7] <- t2-t1
  }
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
    crit4 <- quantile(stat.p4, c(1-alpha))
    TMT <- ifelse(stat4 > crit4, 1, 0)
    return(TMT)
  }
  t1 <- Sys.time()
  tmt.inv <- 0#try(sta4(z, Mp, n, d, alpha), silent=TRUE)
  t2 <- Sys.time()
  if ('try-error' %in% class(tmt.inv)) {
    TMT <- 8888 
    time_m[8] <- 8888
  } else{ 
    TMT <- tmt.inv
    time_m[8] <- t2-t1
  }
  Un_m[8] <- attr(t2-t1,"unit")
  
  ####mdm
  t1 <- Sys.time()
  m.inv <- try(mdm_test(z, c(p,q))$pval, silent=TRUE)  
  t2 <- Sys.time()
  if ('try-error' %in% class(m.inv)) {
    dmdmpVal <- 8888 
    time_m[9] <- 8888
  } else{ 
    dmdmpVal <- ifelse(is.nan(m.inv)|is.na(m.inv), 9999, m.inv)
    time_m[9] <- t2-t1
  }
  Un_m[9] <- attr(t2-t1,"unit")
  
  ######### Han Hallin
  t1 <- Sys.time()
  had.inv <- try({ 
    H_dcovt <- discov.Hallin.stat(x,y)
    for (j in 1:BstrpTimes){
      O1=order(runif(n,0,1))
      YBstrp1 = y[O1,]
      Bs_han[j]  = discov.Hallin.stat(x,YBstrp1)
    } 
  } , silent=TRUE)
  t2 <- Sys.time()
  if ('try-error' %in% class(had.inv)) {
    HaVal <- 8888 
    time_m[10] <- 8888
  }else{ 
    HaVal <- mean(Bs_han >= H_dcovt)
    time_m[10] <- t2-t1
  }
  Un_m[10] <- attr(t2-t1,"unit")
  
  ######Mrcov
  t1 <- Sys.time()
  mrn.inv <- try( { 
    Mrcov <- computestatisticrdcov(x, y)
    for (j in 1:BstrpTimes){
      O2=order(runif(n,0,1));
      YBstrp2 = y[O2,];
      Bs_Mr[j]  = computestatisticrdcov(x, YBstrp2);
    }
  } , silent=TRUE)
  t2 <- Sys.time()
  if ('try-error' %in% class(mrn.inv)) {
    MrVal <- 8888 
    time_m[11] <- 8888
  }else{ 
    MrVal <- mean(Bs_Mr >= Mrcov)
    time_m[11] <- t2-t1
  }
  Un_m[11] <- attr(t2-t1,"unit")
  
  resn <- c(pVal, drpVal, dcortpVal, dhsicpVal,
            Jd, TMT,  dmdmpVal, HaVal, MrVal, sA, sB)
  result <- list(resn = resn, time = time_m, Un_m = Un_m)
  return(result)
}

#####tiems
rep.num <- 2000
mycores <- 40

start <- Sys.time()

res_my <- function(n, p, q, dep, alpha, rep.num, mycores){
res_gauss <- mclapply(1:rep.num, function(rep.sim){
  set.seed(rep.sim+100)
  res_result <- Gauss_T(n, p, q, dep, alpha)
  print(rep.sim)
  write.csv(res_result, paste0("para_res", "_rep.sim", rep.sim, ".csv"))
  return(res_result)
}
, mc.cores = mycores)
end <- Sys.time()
print(end - start)

RE <- matrix(0, 35, rep.num)
for (i in 1:rep.num) {
  RE[,i] <- res_gauss[[i]]$resn
}

Rtime <- matrix(0, 35, rep.num)
for (i in 1:rep.num) {
  Rtime[,i] <- res_gauss[[i]]$time
}


mean_fun <- function(x){
  a <- which((x==8888))
  l <- length(a)
  if (l > 0){
    rm <- mean(x[-a])
  }else{ rm <- mean(x) }
  return(rm)
}
restime <- apply(Rtime, 1, mean_fun)

res1 <- rep(0, 35)
res2 <- rep(0, 35)
fun1 <- function(x){
  a <- length(which((x==8888)))
  b <- length(which((x==9999)))
  x[which((x==8888))] = 1
  x[which((x==9999))] = 1
  y = sum(x[drop=FALSE] <= alpha)/(length(x)-a-b);
  return(y)
}
res1[1:4] <- apply(RE[1:4,], 1, fun1)  

fun2 <- function(x){
  a <- length(which((x==8888)))
  b <- length(which((x==9999)))
  return(a+b)
}
fun3 <- function(x){
  l <- length(which((x==8888)))
  x[which((x==8888))] = 0
  y <- sum(x)/(length(x)- l)
  return(y)
}
res2 <- apply(RE, 1, fun2)

res1[5:8]  <- apply(RE[5:8, ], 1, fun3)
res1[9:11] <- apply(RE[9:11, ], 1, fun1)
res1[12:35] <- apply(RE[12:35, ], 1, mean)

result <- data.frame(list(method=c("pcov", "drcov", "dcor", "hsic", 
                                   "Jdcov", "Jdcov_r", "Jdcov_s", 
                                   "TMT", "mdm", "Hallin", "Mrank",
                                   "gauss-1",   "gauss-3",   "gauss-5", 
                                   "99gauss-1", "99gauss-3", "99gauss-5",
                                   "98gauss-1", "98gauss-3", "98gauss-5", 
                                   "95gauss-1", "95gauss-3", "95gauss-5",
                                   "mgauss-1",   "mgauss-3",   "mgauss-5", 
                                   "m99gauss-1", "m99gauss-3", "m99gauss-5",
                                   "m98gauss-1", "m98gauss-3", "m98gauss-5", 
                                   "m95gauss-1", "m95gauss-3", "m95gauss-5"),
                          prob= res1, Na_nb= res2, times=restime, unit_m =res_gauss[[1]]$Un_m))

return(result)
}



n <- 50
#n_tilde <- floor(n^{6/7})
cat("...........n_", n, "...........\n")

dep <- 0
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, q, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("ex2", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

end <- Sys.time()
print(end - start)


n <- 50
#n_tilde <- floor(n^{6/7})
cat("...........n_", n, "...........\n")
dep <- 1
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, q, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("ex2", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

n <- 50
#n_tilde <- floor(n^{6/7})
cat("...........n_", n, "...........\n")
dep <- 2
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, q, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("ex2", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

n <- 100
#n_tilde <- floor(n^{6/7})
cat("...........n_", n, "...........\n")

dep <- 0
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, q, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("ex2", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

end <- Sys.time()
print(end - start)


n <- 100
#n_tilde <- floor(n^{6/7})
cat("...........n_", n, "...........\n")
dep <- 1
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, q, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("ex2", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

n <- 100
#n_tilde <- floor(n^{6/7})
cat("...........n_", n, "...........\n")
dep <- 2
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, q, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("ex2", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

