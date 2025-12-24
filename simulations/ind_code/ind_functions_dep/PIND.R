library(RcppArmadillo)


# Update this path according to where the source file is located on your system.
Rcpp::sourceCpp('./ind_functions_dep/cv_hat_max_arma.cpp')

#----------------------------------------------------- The proposed independence test -------------------------------------------------------------------------

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




# The proposed independence test with multiplier bootstrap
# Input:
#   x, y:   data matrices with the same number of rows (sample size n)
#   alpha:  significance level in (0, 1)
#   seed:   RNG seed (optional)
#   option: "Rademacher" | "Gaussian" | "Mammen" | "all"
#     "Rademacher" uses Rademacher multipliers;
#     "Gaussian" uses Gaussian multipliers;
#     "Mammen" uses Mammen's multipliers;
#     "all" computes results for all three types of multipliers
#   N: number of bootstrap replications, default 5000

# Output:
#   option:   the chosen option
#   reject:   1: reject the null hypothesis; 0: do not reject the null hypothesis
#   p_value:  p value
#   test_sta: test statistic
#   cv:       critical value
Ind_Gtest <- function(x, y, alpha, seed = 1,
                      option = c("Rademacher", "Gaussian", "Mammen", "all"),
                      N = 5000) {
  option <- match.arg(option)
  
  # ---- ensure x has the larger column dimension ----
  if (ncol(x) < ncol(y)) {
    tmp <- x; x <- y; y <- tmp
  }
  
  # ---- coordinatewise Gaussianization ----
  x_norm <- data_transform(x)$x_trs_my
  y_norm <- data_transform(y)$x_trs_my
  
  # ---- Generate multiplier matrices: N-by-nrow(x) ----
  make_M <- function(kind) {
    B <- N * nrow(x)
    if (kind == "Rademacher") {
      set.seed(seed)
      Mvec <- sample(c(-1, 1), size = B, replace = TRUE, prob = c(0.5, 0.5))  # generate Rademacher multipliers
    } else if (kind == "Gaussian") {
      set.seed(seed)
      Mvec <- rnorm(B, mean = 0, sd = 1)          # generate Gaussian multipliers
    } else if (kind == "Mammen") {
      set.seed(seed)
      Mvec <- rmammen(B, "two-point mass")        # generate Mammen's multipliers
    } else {
      stop("Unknown kind: ", kind)
    }
    matrix(Mvec, nrow = N)       # an N-by-n matrix
  }
  
  # ---- call function cv_arma (cv_hat_max_arma.cpp) for a given multiplier ----
  run_once <- function(kind) {
    M   <- make_M(kind)                  # generate multiplier
    res <- cv_arma(x_norm, y_norm, M, alpha)   
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