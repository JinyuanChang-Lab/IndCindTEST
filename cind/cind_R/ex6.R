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

#Rcpp::sourceCpp('D:/0-xnc-du-dr/1code/0307/cpp/cv_hat_max_arma.cpp')
Rcpp::sourceCpp('cv_hat_max_arma.cpp')

#source('D:/0-xnc-du-dr/CIN/all-functions-CDI.R')
source('all-functions-CDI.R')

#######

Data1 <- function(n, p, m, dep) {
  # 创建学生t分布实例
  q=p
  x <- matrix(rt(n * p, df = 2), nrow = n, ncol = p)
  y <- matrix(rt(n * q, df = 2), nrow = n, ncol = q)
  z <- matrix(rt(n * m, df = 2), nrow = n, ncol = m)
  w <- matrix(rt(n * p, df = 2), nrow = n, ncol = p)
  v <- matrix(rt(n * q, df = 2), nrow = n, ncol = q)
  
  # 获取所有可能的两两列组合的索引
  combinations <- combn(m, 2)
  
  # 计算组合的数量
  num_combinations <- ncol(combinations)
  
  # 对每一对组合的对应元素进行相乘，并将结果存储在x矩阵中
  for (i in 1:num_combinations) {
    col1 <- combinations[1, i]
    col2 <- combinations[2, i]
    x[, i] <- z[, col1] * z[, col2]
    # x[, (i + num_combinations)] <- sin(z[, col1]) + cos(z[, col2]) + z[, col1]^2 + z[, col2]^2
    y[, i] <- z[, col1] + z[, col2]
  }
  
  if (dep != 0) {
    for (i in 1:dep) {
      eps <- rt(n, df = 1)
      x[, i] <- x[, i] + eps + 3 * eps^3
      y[, i] <- y[, i] + eps + 3 * eps^3
    }
  }
  
  # 将数据封装为列表
  data_label <- list(
    x = x,
    y = y,
    z = z
  )
  
  return(data_label)
}

Gauss_CI <- function(n, p, k, alpha){
  dat <- Data1(n,p,k, dep)
  z <- dat$z
  x <- dat$x
  y <- dat$y
  
  sA <- sB <- rep(0, 9)
  time_m <- rep(0, 23)
  Un_m <- rep(0, 23)
  ##################
  
  
  ###########zhu PCD
  t1 <- Sys.time()
  pdcv_inv <- try(pro_Val(x, y, z, alpha),  silent=TRUE)
  #pdcv_inv <- 0
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
  ritv <-  ifelse(rcitp < alpha, 1, 0)
  t2 <- Sys.time()
  time_m[2] <- t2-t1
  Un_m[2] <- attr(t2-t1,"unit")
  ######RcoT
  t1 <- Sys.time()
  rcotp <- RCIT::RCoT(x, y, z)$p
  Rcotv <- ifelse(rcotp < alpha, 1, 0)
  t2 <- Sys.time()
  time_m[3] <- t2-t1
  Un_m[3] <- attr(t2-t1,"unit")
  
  ######
  ### cdsis
  t1 <- Sys.time()
  cdorPval <- cdcsis::cdcov.test(x, y, z)$p.value
  cdor <-  ifelse(cdorPval < alpha, 1, 0)
  t2 <- Sys.time()
  time_m[4] <- t2-t1
  Un_m[4] <- attr(t2-t1,"unit")
  ###1: reject    generalised
  t1 <- Sys.time()
  gcm <- as.numeric(GeneralisedCovarianceMeasure::gcm.test(x, y, z , alpha = 0.05, regr.method="kernel.ridge")$reject)
  t2 <- Sys.time()
  time_m[5] <- t2-t1
  Un_m[5] <- attr(t2-t1,"unit")
  sA[is.na(sA)] <- sB[is.na(sB)] <- 0
  resn <- c(PDCv, ritv, Rcotv, cdor, gcm, sA, sB)
  result <- list(resn = resn, time = time_m, Un_m = Un_m)
}



start <- Sys.time()
res_my <- function(n, p, k, dep, alpha, rep.num, mycores){
  res_gauss <- mclapply(1:rep.num, function(rep.sim){
    res_result <- Gauss_CI(n, p, k, alpha)
    if(rep.sim %% 200 == 0) print(rep.sim)
    #write.csv(res_result, paste0("para_res", "_rep.sim", rep.sim, ".csv"))
    return(res_result)
  }
  , mc.cores = mycores)
  
  RE <- matrix(0, 23, rep.num)
  for (i in 1:rep.num) {
    RE[,i] <- res_gauss[[i]]$resn
  }
  
  Rtime <- matrix(0, 23, rep.num)
  for (i in 1:rep.num) {
    Rtime[,i] <- res_gauss[[i]]$time
  }
  
  
  RE <- matrix(0, 23, rep.num)
  for (i in 1:rep.num) {
    RE[,i] <- res_gauss[[i]]$resn
  }
  
  Rtime <- matrix(0, 23, rep.num)
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
  
  res1 <- rep(0, 23)
  res2 <- rep(0, 23)
  
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
  res1 <- apply(RE, 1, fun3)
  
  
  result <- data.frame(list(method=c("PCD", "PCIT", "PCOT", "CDSIS", "GCD", 
                                     "gauss-1",   "gauss-3",   "gauss-5", 
                                     "99gauss-1", "99gauss-3", "99gauss-5",
                                     "98gauss-1", "98gauss-3", "98gauss-5", 
                                     "mgauss-1",   "mgauss-3",   "mgauss-5", 
                                     "m99gauss-1", "m99gauss-3", "m99gauss-5",
                                     "m98gauss-1", "m98gauss-3", "m98gauss-5"),
                            prob= res1, Nam=res2, times=restime, unit_m =res_gauss[[1]]$Un_m))
return(result)
}



rep.num <- 2000
mycores <- 200
set.seed(67890)



p <- 100
cat("...........p_", p, "...........\n")
k <- 5
alpha <- 0.05


n <- 100
cat("...........n_", n, "...........\n")
dep <- 0
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

end <- Sys.time()
print(end - start)


cat("...........n_", n, "...........\n")
dep <- p/10
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/5
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

n <- 200
cat("...........n_", n, "...........\n")
dep <- 0
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/10
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/5
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))


p <- 400
cat("...........p_", p, "...........\n")
k <- 5
alpha <- 0.05

n <- 100
cat("...........n_", n, "...........\n")
dep <- 0
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/10
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/5
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

n <- 200
cat("...........n_", n, "...........\n")
dep <- 0
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/10
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/5
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

p <- 1600
cat("...........p_", p, "...........\n")
k <- 5
alpha <- 0.05

n <- 100
cat("...........n_", n, "...........\n")
dep <- 0
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/10
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/5
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

n <- 200
cat("...........n_", n, "...........\n")
dep <- 0
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

end <- Sys.time()
print(end - start)


cat("...........n_", n, "...........\n")
dep <- p/10
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))

cat("...........n_", n, "...........\n")
dep <- p/5
cat("...........dep_", dep, "...........\n")
#############dep=0##########
result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
write.csv(result_m, paste0("cex6", "_n", n, "_p", p, "_a", alpha*100,"_dep",dep, ".csv"))
