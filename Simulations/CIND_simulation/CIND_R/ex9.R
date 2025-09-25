rm(list=ls()) 
##########  conditional independece test   ###############
library(MASS)
library(parallel)
library(glmnet)
library(cdcsis)
library(GeneralisedCovarianceMeasure)
#library(devtools)
library(remotes)
#remotes::install_github("ericstrobl/RCIT", force = TRUE)
library(RCIT)

source('cind-functions.R') 
source('PCIND.R')

Rcpp::sourceCpp('cv_hat_max_arma.cpp')
 


Data1 <- function(n, p, k, dep){
  q=p
  x <- matrix(rnorm(n*p), n, p)
  y <- matrix(rnorm(n*p), n, q)
  z <- matrix(rnorm(n*k), n, k)
  w <- matrix(rnorm(n*k), n, k)
  v <- matrix(rnorm(n*k), n, k)
  for (i in 1:k) {
    a <- 0.5*(z[, i]^3/7 + z[, i]/2)
    b <- (z[, i]^3/2 + z[, i])/3
    a1 <- a + tanh(w[ ,i])
    x[, i] <- a1 + a1^3/3
    b1 <- b + v[ ,i]
    y[, i] <- b1 +  tanh(b1/3)
  }
  if(dep!=0){
    for (i in 1:dep) {
      epsilon <- rnorm(n, 0, 1)
      x[, i] <- 0.5 * x[, i] + 3*cosh(epsilon)
      y[, i] <- 0.5 * y[, i] + 3*cosh(epsilon^2)
    }
  }
  return(list(x=x,y=y,z=z))
}




start <- Sys.time()

# Wrapper function to run Gauss_CI function (in all-functions-CDI.R) in parallel and aggregate the results
res_my <- function(n, p, k, dep, alpha, rep.num, mycores){ 
  # rep.num: the number of replications
  # mycores: the number of cores used for parallelization
  
  res_gauss <- mclapply(1:rep.num, function(rep.sim){                     # Use mclapply for parallel computation 
    if(rep.sim %% 200 == 0)  cat('---- Test times = ', rep.sim, '----\r',file = "ex9log.txt", append=TRUE)
    set.seed(rep.sim+100)                                                 # Set a replication-specific seed for reproducibility
    res_result <- Gauss_CI(n, p, k, alpha, rep.sim)                                  
    if(rep.sim %% 200 == 0) print(rep.sim)                                  
    #write.csv(res_result, paste0("para_res", "_rep.sim", rep.sim, ".csv"))  
    return(res_result)                                                     
  }
  , mc.cores = mycores, mc.preschedule = FALSE)                                                     
  #On Windows system, the number of cores should be set to 1. 
  #On macOS, Linux, or other Unix-like systems, the number of cores can be adjusted as appropriate. 
  
  #results of the proposed conditional independence test based on linear regressions and the competing methods
  RE <- matrix(0, 23, rep.num)                                              
  for (i in 1:rep.num) {                                                    
    RE[,i] <- res_gauss[[i]]$resn                                          
  }
  
  ###### running time of all methods 
  Rtime <- matrix(0, 23, rep.num)                                           
  for (i in 1:rep.num) {                                                    
    Rtime[,i] <- res_gauss[[i]]$time                                        
  }
  
  # Compute the proportion of rejections of the null across the valid replications
  mean_fun <- function(x){                                                 # x is a rep.num-by-1 vector with entries in {0,1}; 1 indicates rejection of the null hypothesis  
    a <- which((x==8888))                                                  # find indices of invalid runs
    l <- length(a)                                                         # count the number of invalid runs
    if (l > 0){                                                            # take the average excluding invalid results
      rm <- mean(x[-a])                                                     
    }else{ rm <- mean(x) }                                                  
    return(rm)                                                              
  }
  
  restime <- apply(Rtime, 1, mean_fun)                                     # average running time of each method 
  
  res1 <- rep(0, 23)                                                        
  res2 <- rep(0, 23)                                                        
  
  # count the number of invalid runs
  fun2 <- function(x){                                                      
    a <- length(which((x==8888)))                                          
    b <- length(which((x==9999)))                                           
    return(a+b)                                                             
  }
  
  # compute the proportion of rejections of the null across the valid replications
  fun3 <- function(x){                                                     # x is a rep.num-by-1 vector with entries in {0,1}; 1 indicates rejection of the null hypothesis  
    l <- length(which((x==8888)))                                          # count the number of invalid runs
    x[which((x==8888))] = 0                                                   
    y <- sum(x)/(length(x)- l)                                             # Average over valid entries
    return(y)                                                               
  }
  ####count the number of invalid runs
  res2 <- apply(RE, 1, fun2) 
  
  #### calculate the proportion of rejections of the null for each methods
  res1 <- apply(RE, 1, fun3)                                                
  
  
  res1_fin = c(res1[6:8], res1[5],res1[1],res1[2:3],res1[4])
  res2_fin = c(res2[6:8], res2[5],res2[1],res2[2:3],res2[4])
  result <- data.frame(list(method=c("gaussian","mammen","rademacher", "GCM", "PCD", "PCIT", "PCOT", "cdCov"),
                            prob= res1_fin, Nam=res2_fin))
  
  cat("---- result=\n", file = "ex9log.txt", append = TRUE)
  capture.output(print(result), file = "ex9log.txt", append = TRUE)
  cat("\n----\n", file = "ex9log.txt", append = TRUE)
  return(result)
}



rep.num <- 2000
mycores <- 100

k <- 5
alpha <- 0.05

p <- 100
cat("...........p_", p, "...........\n")

for (n in c(100)) {
  for (dep in c(0, p/10, p/5)) {
    cat("...........n_", n, "...........\n")
    cat("...........dep_", dep, "...........\n")
    
    start <- Sys.time()
    result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
    
    write.csv(
      result_m,
      paste0("cex9", "_n", n, "_p", p, "_a", alpha*100, "_dep", dep, ".csv"),
      row.names = FALSE
    )
    
    end <- Sys.time()
    print(end - start)
  }
}