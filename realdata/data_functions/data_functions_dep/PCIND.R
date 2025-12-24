library(RcppArmadillo)
library(glmnet)

Rcpp::sourceCpp('./data_functions_dep/cv_hat_max_arma.cpp')

#--------------------------------------- The proposed conditional independence test based on linear regressions ------------


######################    Coordinatewise   Gaussianization  ###############
#   - This is a column-wise, marginal Gaussianization.
data_transform <- function(x){
  n <- nrow(x)                            # sample size
  p <- ncol(x)                            # number of variables
  x_t2 <- matrix(0, nrow = n, ncol = p)   # generate n-by-p matrix
  
  for (j in 1:p) {                   # process each column independently
    # Calculate the empirical distribution function, scaled by n/(n+1) 
    ecdf_values <- ecdf(x[, j])(x[, j]) * n / (n + 1)
    x_t2[, j] <- ecdf_values
    
    for (i in 1:n) {                      
      if (x_t2[i, j] != 0 && x_t2[i, j] != 1) {
        x_t2[i, j] <- qnorm(x_t2[i, j])  # standard normal quantile at x_t2[i, j]
      } else {
        x_t2[i, j] <- 0                  # avoid Inf values 
      }
    }
  }
  
  list(x_trs = x, x_trs_my = x_t2)    
}



################# generate Mammen's multiplier: construct = "two-point mass" ################
rmammen <- function(n,
                    construct = c("normal-2", "normal-1", "two-point mass")){
  construct <- match.arg(construct)
  
  if (length(n) > 1) n <- length(n)
  
  if (construct == "two-point mass"){
    .vals <- c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2) 
    .probs <- rev(abs(.vals)/sqrt(5))
    
    sample(.vals, size = n, replace = TRUE, prob = .probs)    #generate random samples taking values ±(sqrt(5)±1)/2 with probabilities ensuring mean 0 and variance 1.
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



######## penalty for residual (Lasso-based penalty) ####################


fun_penalty <- function(x, z){                        # Get residuals by fitting Lasso of x (1-dimensional) on z (p-dimensional)
  set.seed(12345)
  p <- dim(z)[2]                                      # number of predictors (columns of z)
  obj.cv1 <- cv.glmnet(z,x, parallel = FALSE)        
  lambda_my <-obj.cv1$lambda.min                      # select lambda that minimizes the CV error
  
  
  beta_my <- glmnet(z, x, lambda = lambda_my)$beta    # Fit Lasso at chosen lambda and extract sparse coefficients 
  index_nz <- beta_my@i                               # row indices of nonzero coefficients in sparse matrix (zero-based indexing)
  num <- beta_my@x                                    # corresponding nonzero coefficient values
  beta <- rep(0, p)
  beta[index_nz+1] <- num
  resid  <- x - z %*% beta                            # calculate the residuals
  return(resid)
}



# The proposed conditional independence test based on linear regressions with multiplier bootstrap 
# Input:
#   x, y, z: data matrices with the same number of rows (sample size n)
#   alpha:   significance level in (0, 1)
#   seed:    RNG seed (optional)
#   option:  "Rademacher" | "Gaussian" | "Mammen" | "all"
#     "Rademacher" uses Rademacher multipliers;
#     "Gaussian" uses Gaussian multipliers;
#     "Mammen" uses Mammen's multipliers;
#      "all" computes results for all three types of multipliers
#   N:    number of bootstrap replications, default 5000

# Output:
#   option:   the chosen option
#   reject:   1: reject the null hypothesis; 0: do not reject the null hypothesis
#   p_value:  p value
#   test_sta: test statistic
#   cv:       critical value
Cind_Gtest <- function(x, y, z, alpha, seed = 1,
                      option = c("Rademacher", "Gaussian", "Mammen", "all"),
                      N = 5000) {
  option <- match.arg(option)
  
  # ---- ensure x has the larger column dimension ----
  if (ncol(x) < ncol(y)) {
    tmp <- x; x <- y; y <- tmp
  }
  
  # ---- coordinatewise Gaussianization ----
  z_norm <- data_transform(z)$x_trs_my                                     
  x_norm <- data_transform(x)$x_trs_my                                     
  y_norm <- data_transform(y)$x_trs_my                                     
  if(dim(z)[2]==1){                                                       # if z is univariate, use simple linear regression
    Epsilon <- matrix(lm(x_norm~z_norm-1)$residuals, n)                   # residuals of x_norm on z_nrom  
    Delta <- matrix(lm(y_norm~z_norm-1)$residuals, n)                     # residuals of y_norm on z_norm  
  }else{                                                                   
    Epsilon <- apply(x_norm, 2, fun_penalty, z_norm)                      # apply Lasso regression of each column of x_norm on z_norm, and get residuals 
    Delta   <- apply(y_norm, 2, fun_penalty, z_norm)                      # apply Lasso regression of each column of y_norm on z_norm, and get residuals 
  }
  
  # ---- Generate multiplier matrices: N-by-nrow(x) ----
  make_M <- function(kind) {
    B <- N * nrow(x)
    if (kind == "Rademacher") {
      set.seed(seed)
      Mvec <- sample(c(-1, 1), size = B, replace = TRUE, prob = c(0.5, 0.5))
    } else if (kind == "Gaussian") {
      set.seed(seed)
      Mvec <- rnorm(B, mean = 0, sd = 1)
    } else if (kind == "Mammen") {
      set.seed(seed)
      Mvec <- rmammen(B, "two-point mass")
    } else {
      stop("Unknown kind: ", kind)
    }
    matrix(Mvec, nrow = N)   
  }
  
  # ---- call function cv_arma (cv_hat_max_arma.cpp) for a given multiplier ----
  run_once <- function(kind) {
    M   <- make_M(kind)                  # generate multiplier
    res <- cv_arma(Epsilon, Delta, M, alpha)   
    tsta <- res$ts                       # test statistic
    cv   <- res$cv_e                     # critical value
    reject <- as.numeric(tsta > cv)      # reject the null hypothesis (=1) if the test statistic > critical value
    list(
      option   = kind,
      reject   = reject, 
      p_value  = res$p_value,            # p value
      test_sta = tsta,
      cv       = cv
    )
  }
  
  if (option == "all") {
    out_r <- run_once("Rademacher")
    out_g <- run_once("Gaussian")
    out_m <- run_once("Mammen")
    return(list(
      Rademacher = out_r,
      Gaussian   = out_g,
      Mammen     = out_m
    ))
  } else {
    return(run_once(option))
  }
}