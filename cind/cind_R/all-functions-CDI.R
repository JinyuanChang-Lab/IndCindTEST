#########functions of CDI
library(glmnet)
data_transform <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  x_t2 <- matrix(0, nrow = n, ncol = p)
  
  for (j in 1:p) {
    ecdf_values <- ecdf(x[, j])(x[, j]) * n / (n + 1)
    x_t2[, j] <- ecdf_values
    
    for (i in 1:n) {
      if (x_t2[i, j] != 0 && x_t2[i, j] != 1) {
        x_t2[i, j] <- qnorm(x_t2[i, j])
      } else {
        x_t2[i, j] <- 0
      }
    }
  }
  
  list(x_trs = x, x_trs_my = x_t2)
}

########penalty for residual
fun_penalty_para <- function(x, z, index){
  p <- dim(z)[2]
  obj.cv1 <- cv.glmnet(z,x[, index])
  lambda_my <-obj.cv1$lambda.min 
  beta_my <- glmnet(z, x[, index], lambda = lambda_my)$beta
  ind_nz <- beta_my@i
  num <- beta_my@x
  beta <- rep(0, p)
  beta[ind_nz+1] <- num
  resid  <- x[, index] - z %*% beta
  resid <- matrix(resid, n)
  return(resid)
}

fun_penalty <- function(x, z){
  p <- dim(z)[2]
  obj.cv1 <- cv.glmnet(z,x, parallel = FALSE)
  lambda_my <-obj.cv1$lambda.min 
  beta_my <- glmnet(z, x, lambda = lambda_my)$beta
  index_nz <- beta_my@i
  num <- beta_my@x
  beta <- rep(0, p)
  beta[index_nz+1] <- num
  resid  <- x - z %*% beta
  return(resid)
}

########zhu PCD

CDFNonParamv<-function(X,Y,X0,Y0,bw,family = c("Normal", "Epan", "BOX")){
  n <- dim(X)[1];
  p <- dim(X)[2];
  m <- dim(Y)[1];
  d <- dim(Y)[2];
  
  n0 <- dim(X0)[1];
  p0 <- dim(X0)[2];
  m0 <- dim(Y0)[1];
  d0 <- dim(Y0)[2];
  if (n0 == 1 && p0 == 1)
  {X0 <- rep(X0,p)}
  
  # Default window parameter is optimal for normal distribution
  if (is.null(bw)){
    sig <- apply(X,2,mad)
    bw  <- sig * (4/(3*n))^(1/(p+4))
  }
  bw<-as.vector(bw)
  ell <- length(bw)
  if (ell==1){
    bw<-as.vector(rep(bw,p))
  }
  if (length(bw) != p) 
    stop("Length of h must equal to the number of predictors")
  
  KerUni<-array(0,dim=c(n0,n,p))
  KerWgt=matrix(1,n0,n)
  for (dim in 1:p){
    Xt <- X[,dim]
    Xt0 <- X0[,dim]
    U<-(matrix(rep(Xt0,n),n,n)-t(matrix(rep(Xt,n0),n0,n0)))/bw[dim]
    if (family == "Normal" ){
      KerUni[,,dim]= exp(-0.5 * U^2) / sqrt(2*pi);
    }else if (family == "Epan" ){
      KerUni[,,dim]= 0.75 * (1 - U^2) * (abs(U) <= 1);
    }else if (family == "BOX" ){
      KerUni[,,dim]= (abs(U)<=a) / (2 * a);
    }
    KerWgt= KerUni[,,dim]*KerWgt
  }
  KerSum<-apply(KerWgt,2,sum)
  Ywgt<-matrix(1,m0,m)
  for (dim in 1:d){
    Yt <- Y[,dim]
    Yt0<- Y0[,dim]
    Id <- matrix(rep(Yt,m0),m0,m0)<=t(matrix(rep(Yt0,m),m,m))
    Ywgt<-Id*Ywgt}
  CDF=KerWgt%*%Ywgt/matrix(rep(KerSum,m0),m0,n)
  return(CDF)  
}

ProConD<-function(x,y,z,zn,hh){
  # Measure Dependence between a vector Y and a vector X conditional on Z 
  n <- dim(y)[1];
  d <- dim(y)[2];
  indy<-matrix(1,n,n)
  for (dim in 1:d){
    Y0 <- y[,dim]
    Id <- matrix(rep(Y0,n),n,n)<=t(matrix(rep(Y0,n),n,n))
    indy<-Id*indy}
  Fyoxb<- CDFNonParamv(zn,y,zn,y,hh,family ='Normal')
  diff<-indy-Fyoxb
  cov0<-kronecker(matrix(1,n,1),diff)*matrix(rep(diff,each=n),n*n,n)
  mcov<-apply(cov0, 1, mean) 
  b<-matrix(mcov, nrow = n, ncol = n)
  bjbar<-matrix(rep(apply(b,2,mean),each=n),n,n)
  bbar<-mean(apply(b,2,mean))
  bb<-b-bjbar-t(bjbar)+bbar
  
  xxup<-cbind(x,z)%*%t(cbind(x,z))+1
  xxdown<-sqrt(diag(xxup))
  xxdown<-xxdown%*%t(xxdown)
  
  rhok<-1/4+ asin(pmin(pmax((xxup/xxdown),-1.0), 1.0))/(2*pi)
  rhok<-abs(rhok)
  
  bbrho<-bb*rhok
  conCov<-mean(apply(bbrho,2,mean))
  evar<-mean(apply(diff,2,var)*(n-1)/n)
  conCor <-sqrt(conCov/evar)
  return(conCor)
  
}

pro_Val <- function(x, y, z, alpha){
  n <- dim(x)[1]
  h0 <-  n^(-1/3)
  h <- 1.06 * h0 * sd(x)
  pSta <- ProConD(x,y,z,z,h)
  Bstrp <- 500
  STb <- rep(0, Bstrp)
  for (t in 1:Bstrp) {
    index <- sample(seq(1, n, 1), n, replace= TRUE, prob = rep(1/n, n))                       
    xtstar <- matrix(0, dim(x)[1], dim(x)[2])
    for (i in 1: dim(x)[2]) {
      xtstar[, i] <- x[index, i] + h * rnorm(n, 0,1) 
    }
    hstar <- 1.06 * h0 * sd(xtstar)
    k02 <- array(0 ,dim = c(n,n,dim(x)[2]))
    for (j in 1:dim(x)[2]) {
      k02[, , j] <- (x[, j] %*% matrix(1, 1, n) - matrix(1, n, 1) %*% t(x[, j]) )/h
    } 
    k02 <- (2*pi)^(-1/2)* exp(- k02^2/2)/h
    m <- k02[,,1]
    
    k02_new <- array(k02, dim=c(n*n, 1, dim(x)[2]))
    a <- matrix(apply(k02_new, 1, prod), n, n)
    b <- apply(a, 2, sum)
    
    w0 <- a / (matrix(1, n, 1) %*% matrix(b,1))
    s <- apply(w0, 2, cumsum)
    ss1 <- rbind(s, runif(n,0,1))
    ss2 <- rbind(s, runif(n,0,1))
    fun_max <- function(x){
      x_new <- order(x)
      d <- which(x_new== length(x))
      return(d)
    }
    i1 <- apply(ss1, 2, fun_max)
    i2 <- apply(ss2, 2, fun_max)
    ytstar=matrix(y[i1,],n)
    ztstar=matrix(z[i2,],n)    
    pSta_b <- ProConD(xtstar,ytstar,ztstar,ztstar,hstar)
    STb[t] <- pSta_b
    #print(t)
  }
  proVal <- ifelse(pSta > quantile(STb, 1-alpha), 1, 0)
  return(proVal)
}
