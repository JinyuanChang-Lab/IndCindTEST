rm(list=ls()) 
##########  independece test   ###############
 
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


source('ind-functions.R')
source('PIND.R')

Rcpp::sourceCpp('cv_hat_max_arma.cpp')
 



#####generate data-Example 3 
Data1 <- function(n, p, q, dep){
  if(dep==0){
    x <- matrix(runif(n*p, 0, 2*pi), n, p)
    y <- matrix(runif(n*q, 0, 2*pi), n, q)
  }
  else{
    if(dep==1){d <- p/20
    }else if(dep==2){d <- p/10
    }else{d <- p/5}
    w <- matrix(runif(n*d, 0, 2*pi), n, d)
    x1 <- (sin(w))^2 
    y1 <- (cos(w))^2
    x <- cbind(x1, matrix(runif(n*(p-d), 0, 2*pi), n))
    y <- cbind(y1, matrix(runif(n*(q-d), 0, 2*pi), n))
  }
  dat <- list(x=x, y=y)
  return(dat)  
}




start <- Sys.time()

# Wrapper function to run Gauss_T function (in all-functions.R) in parallel and aggregate the results
res_my <- function(n, p, q, dep, alpha, rep.num, mycores){   
  # rep.num: the number of replications
  # mycores: the number of cores used for parallelization
  
  res_gauss <- mclapply(1:rep.num, function(rep.sim){      # Use mclapply for parallel computation             
    set.seed(rep.sim+100)                                  # Set a replication-specific seed for reproducibility
    
    if(rep.sim %% 200 == 0)  cat('---- Test times = ', rep.sim, '----\r',file = "ex3log.txt", append=TRUE)
    res_result <- Gauss_T(n, p, q, dep, alpha, rep.sim)                             
    print(rep.sim)                                                          
    #write.csv(res_result, paste0("para_res", "_rep.sim", rep.sim, ".csv")) 
    return(res_result)                                                      
  }
  , mc.cores = mycores, mc.preschedule = FALSE)                             
  #On Windows system, the number of cores should be set to 1. 
  #On macOS, Linux, or other Unix-like systems, the number of cores can be adjusted as appropriate. 
  end <- Sys.time()                                                         
  print(end - start)                                                        
  
  #results of the proposed independence test, and competing methods without coordinatewise Gaussianization
  RE <- matrix(0, 34, rep.num)                                              
  for (i in 1:rep.num) {                                                    
    RE[,i] <- res_gauss[[i]]$resn                                          # the results without coordinatewise Gaussianization of the i-th replication 
  }
  
  #results of the proposed independence test, and competing methods with coordinatewise Gaussianization 
  RE_gaussian <- matrix(0, 34, rep.num)                                     
  for (i in 1:rep.num) {                                                    
    RE_gaussian[,i] <- res_gauss[[i]]$resn_gaussian                        # the results with coordinatewise Gaussianization of the i-th replication 
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
  
  
  res1 <- rep(0, 34)                                                        
  res1_gaussian <- rep(0, 34)                                               
  res2 <- rep(0, 34)                                                        
  res2_gaussian <- rep(0, 34)  
  
  # Compute the proportion of rejections of the null across the valid replications
  fun1 <- function(x){                                                     # x is a rep.num-by-1 vector containing p-values
    a <- length(which((x==8888)))                                          # find indices of invalid runs
    b <- length(which((x==9999)))                                          # find indices of invalid runs
    x[which((x==8888))] = 1                                                 
    x[which((x==9999))] = 1                                                 
    y = sum(x[drop=FALSE] <= alpha)/(length(x)-a-b);                       # calculate the proportion of p-values <= alpha among valid runs
    return(y)                                                               
  }
  
  #### calculate the proportion of rejections of the null for Pcor, rdCov, dCor, dHSIC methods
  res1[1:4] <- apply(RE[1:4,], 1, fun1)                                    # without coordinatiwise Gaussianization
  res1_gaussian[1:4] <- apply(RE_gaussian[1:4,], 1, fun1)                  # with coordinatiwise Gaussianization
  
  
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
  res2 <- apply(RE, 1, fun2)                                               # without coordinatiwise Gaussianization
  res2_gaussian <- apply(RE_gaussian, 1, fun2)                             # with coordinatiwise Gaussianization
  #### calculate the proportion of rejections of the null for JdCov_R method (res1[6])
  res1[5:7]  <- apply(RE[5:7, ], 1, fun3)                     # without coordinatiwise Gaussianization
  res1_gaussian[5:7]  <- apply(RE_gaussian[5:7, ], 1, fun3)   # with coordinatiwise Gaussianization
  #### calculate the proportion of rejections of the null for GdCov, Hallan, mrdCov  methods
  res1[8:10] <- apply(RE[8:10, ], 1, fun1)                     # without coordinatiwise Gaussianization
  res1_gaussian[8:10] <- apply(RE_gaussian[8:10, ], 1, fun1)   # with coordinatiwise Gaussianization
  #### calculate the proportion of rejections of the null for the proposed method (res1[11])
  res1[11:34] <- apply(RE[11:34, ], 1, mean)                         # without coordinatiwise Gaussianization        
  res1_gaussian[11:34] <- apply(RE_gaussian[11:34, ], 1, mean)       # with coordinatiwise Gaussianization
  
  ## aggregate the results
  res1_fin <- c(res1[11:13], res1[1:4], res1[6], res1[8:10]) 
  res1_gaussian_fin <- c(res1_gaussian[11:13], res1_gaussian[1:4], res1_gaussian[6], res1_gaussian[8:10]) 
  res2_fin <- c(res2[11:13], res2[1:4], res2[6], res2[8:10]) 
  res2_gaussian_fin <- c(res2_gaussian[11:13], res2_gaussian[1:4], res2_gaussian[6], res2_gaussian[8:10]) 
  result <- data.frame(list(method=c("gaussian","mammen", "rademacher", "Pcor", "rdCov", "dCor", "dHSIC", "JdCovR", "GdCov", "Hallin", "mrdCov"),
                            prob= res1_fin, prob_gaussian= res1_gaussian_fin, Na_nb= res2_fin,Na_nb_gaussian= res2_gaussian_fin))
  #cat('---- result= ', result, '----\r',file = "ex1log.txt", append=TRUE)
  cat("---- result=\n", file = "ex3log.txt", append = TRUE)
  capture.output(print(result), file = "ex3log.txt", append = TRUE)
  cat("\n----\n", file = "ex3log.txt", append = TRUE)
  return(result)
}





#####tiems
rep.num <- 2000
mycores <- 100


p=100 
q=100 
alpha=0.05

nall <- c(50)    # sample sizes (can extend, e.g., c(50, 100, 200))

for (n in nall) {
  
  for (dep in c(0,1, 2)) {
    
    cat("...........n_", n, "...........\n")
    cat("...........dep_", dep, "...........\n")
    
    start <- Sys.time()
    result_m <- res_my(n, p, q, dep, alpha, rep.num, mycores)
    
    
    write.csv(
      result_m,
      paste0("ex3", "_n", n, "_p", p, "_a", alpha*100, "_dep", dep, ".csv"),
      row.names = FALSE
    )
    
    end <- Sys.time()
    print(end - start)
  }
}





