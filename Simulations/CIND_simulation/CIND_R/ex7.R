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



Data1<- function(n, p, m, dep) {
  # Initialize dimensions: set q = p for convenience
  q=p
  
  # Initialize matrices:
  # x, y ~ i.i.d. N(0,1); z ~ i.i.d. Uniform(-1,1)
  
  x <- matrix(rnorm(n * p, 0, 1), nrow = n, ncol = p)
  y <- matrix(rnorm(n * q, 0, 1), nrow = n, ncol = q)
  z <- matrix(runif(n * m) * 2 - 1, nrow = n, ncol = m)  
  
  
  # Convert the integer 'dep' to a correlation-like parameter rho in (-1,1)
  rho <- dep * 0.1
  
  # Compute beta from rho; requires |rho| < 1
  beta <- rho / (2 * sqrt(1 - rho^2))
  
  # Generate the first m columns of y and x using z[:, i]
  for (i in 1:m) {
    # Create mean-zero noise for y by summing 48 uniforms in (-0.25, 0.25)
    eps_y1 <- matrix(runif(n * 48) * 0.5 - 0.25, nrow = n, ncol = 48)
    eps_y <- rowSums(eps_y1)  
    
    # Nonlinear function of z plus noise
    y[, i] <- z[, i] + 0.25 * z[, i]^2 + eps_y
    
    # Gaussian noise for x
    eps_x <- rnorm(n,0,1)  
    # x depends on y and z (linear in y, linear in z) plus Gaussian noise
    x[, i] <- 5 * beta * y[, i] + z[, i] + eps_x
  }
  
  # For remaining columns (i = m+1, ..., p), update x using y[, i]
  # (Here y[, i] exists because y has p columns.)
  for (i in (m + 1):p) {
    x[, i] <- x[, i] + 5 * beta * y[, i] + rnorm(n)
  }
  return(list(x=x,y=y,z=z))
}





start <- Sys.time()

# Wrapper function to run Gauss_CI function (in all-functions-CDI.R) in parallel and aggregate the results
res_my <- function(n, p, k, dep, alpha, rep.num, mycores){ 
  # rep.num: the number of replications
  # mycores: the number of cores used for parallelization
  
  res_gauss <- mclapply(1:rep.num, function(rep.sim){                     # Use mclapply for parallel computation 
    if(rep.sim %% 200 == 0)  cat('---- Test times = ', rep.sim, '----\r',file = "ex7log.txt", append=TRUE)
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
  
  cat("---- result=\n", file = "ex7log.txt", append = TRUE)
  capture.output(print(result), file = "ex7log.txt", append = TRUE)
  cat("\n----\n", file = "ex7log.txt", append = TRUE)
  return(result)
}




rep.num <- 2000
mycores <- 100



k <- 5
alpha <- 0.05

p <- 100
cat("...........p_", p, "...........\n")

for (n in c(100)) {
  for (dep in c(0, 7, 8)) {
    cat("...........n_", n, "...........\n")
    cat("...........dep_", dep, "...........\n")
    
    start <- Sys.time()
    result_m <- res_my(n, p, k, dep, alpha, rep.num, mycores)
    
    write.csv(
      result_m,
      paste0("cex7", "_n", n, "_p", p, "_a", alpha*100, "_dep", dep, ".csv"),
      row.names = FALSE
    )
    
    end <- Sys.time()
    print(end - start)
  }
}