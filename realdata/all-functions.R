###packages
library(clue)
library(adagio)
library(energy)
library(pracma)
library(randtoolbox)



######################       Gaussian  ###############
#data transform
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

######arma.cpp

######################       zhu  pCov    ####################
Subtractr=function(Z,r){
  n=nrow(Z);p=ncol(Z);
  Zctr =Z - matrix(1,n,1)%*%matrix(Z[r,],1,p);                        #  Z(i,:) - Z(r,:)
  Zctr=Zctr[-r,];
  
  
  ZrNorm =matrix(sqrt(rowSums(Zctr^2)),n-1,1);                          #norm|Z(i,:) - Z(r,:)|
  ZrStd =Zctr/(ZrNorm%*%matrix(1,1,p));                                 #Z(i,:) - Z(r,:)/|Z(i,:) - Z(r,:)|
  A=ZrStd%*%t(ZrStd); 
  A[A>1]=1
  A[A<(-1)]=-1
  A=acos(A)
  A[is.nan(A)]=0;
  return(A)                                            
}

My_two_MULTI_U=function(X,Y){
  
  n=nrow(X);p=ncol(X); q=ncol(Y); ##nrow(X)=nrow(Y)
  IN = rep(0,n);
  SA = rep(0,n);
  
  for (r in 1:n) {
    
    MX=Subtractr(X,r);
    MY=Subtractr(Y,r);
    XYn2=MX%*%MY
    
    diagMX<-rep(diag(MX),each=n-2)
    diagMY<-rep(diag(MY),each=n-2)
    offdiagM<-as.logical(lower.tri(MX)+upper.tri(MX))
    VecOffdiagMX<-matrix(MX,(n-1)^2,1)[offdiagM]
    VecOffdiagMY<-matrix(MY,(n-1)^2,1)[offdiagM]
    
    Sr3=(sum(XYn2)-sum(diag(XYn2))-sum(diagMX*VecOffdiagMY)-sum(diagMY*VecOffdiagMX))/((n-1)*(n-2)*(n-3))
    
    IN[r] =(sum(MX*MY)-sum(diag(MX*MY)))/((n-1)*(n-2))+((sum(MX)-sum(diag(MX)))/((n-1)*(n-2)))*((sum(MY)-sum(diag(MY)))/((n-1)*(n-2)))-2*Sr3
    SA[r]=((sum(MX)-sum(diag(MX)))/((n-1)*(n-2)))*((sum(MY)-sum(diag(MY)))/((n-1)*(n-2)));
  }
  #Indep_Index= mean(IN);
  Indep_Index= mean(IN)/(pi^2-mean(SA));
  return(Indep_Index)
}

################ zhang Jdcov   ###############
v.center <- function(x){
  if (is.matrix(x)) {n=dim(x)[1]
  if (isSymmetric(x)) A=x else A=as.matrix(dist(x))} else {n=length(x); A=as.matrix(dist(x))}
  R=rowSums(A)
  C=colSums(A)
  T=sum(A)                             
  r=matrix(rep(R,n),n,n)/n
  c=t(matrix(rep(C,n),n,n))/n
  t=matrix(T/n^2,n,n)
  UA=-(A-r-c+t)
  return(UA)
}

Jdcov.sq.V <- function(x,c){    ## Computes the V-Statistic type estimator of JdCov
  n=dim(x)[1]
  d=dim(x)[2]
  A <- v.center(x[,1]) + c 
  for(i in 2:d) A <- A*( v.center(x[,i]) + c )
  return( sum(A)/n^2 - c^d )
}

u.center <- function(x){
  if (is.matrix(x)) {n=dim(x)[1]
  if (isSymmetric(x)) A=x else A=as.matrix(dist(x))} else {n=length(x); A=as.matrix(dist(x))}
  R=rowSums(A)
  C=colSums(A)
  T=sum(A)
  r=matrix(rep(R,n),n,n)/(n-2)
  c=t(matrix(rep(C,n),n,n))/(n-2)
  t=matrix(T/(n-1)/(n-2),n,n)
  UA=-(A-r-c+t)
  diag(UA)=0
  return(UA)
}

Jdcov.sq.U <- function(x,c){    ## Computes the U-Statistic type estimator of JdCov
  n=dim(x)[1]
  d=dim(x)[2]
  A <- u.center(x[,1]) + c 
  for(i in 2:d) A <- A*( u.center(x[,i]) + c )
  return( sum(A)/n/(n-3)-c^d*n/(n-3) )
}

dCov.sq.U <- function(x,n) {
  A <- u.center(x)
  return( sum(A^2)/n/(n-3) )
}

Jdcov.sq.US <- function(x,c){   ## Computes the bias-corrected estimator of JdCov_S
  n=dim(x)[1]
  d=dim(x)[2]
  A <- u.center(x[,1])/sqrt( dCov.sq.U(x[,1],n) ) + c         
  for(i in 2:d) A <- A*( u.center(x[,i])/sqrt( dCov.sq.U(x[,i],n) ) + c )
  return( sum(A)/n/(n-3)-c^d*n/(n-3) )
}

###################    Han   Hallin ###########
do <- function(p) { x = rnorm(p); return(x/sqrt(sum(x^2))) }

do.grid <- function(n,p){
  if (p == 1) {
    nR = n %/% 2; n0 = n %% 2
    g = matrix(c(-(nR:1)/(nR+1), rep(0,n0), (1:nR)/(nR+1)), n)
  } else {
    nR = floor(n^0.25); nS = n %/% nR; n0 = n %% nR
    g = rbind(((1:nR)/(nR+1)) %x% t(replicate(nS, do(p))), matrix(0,n0,p))
  }
  return(g)
}

Hallin.rank <- function(X, const=1e4, constC=1e8){
  n = nrow(X); p = ncol(X)
  X.Hallin = do.grid(n,p)
  if (p == 1){
    X.rank = X.Hallin[rank(X),]
  } else {
    X.DIST = matrix(NA,n,n)
    for (u in 1:n){
      for (v in 1:n){
        X.DIST[u,v] = sum((X[u,] - X.Hallin[v,])^2)
      }
    }
    if (const*max(X.DIST) < constC) {
      X.rank = X.Hallin[assignment(round(const*X.DIST))$perm,]
    } else {
      X.rank = X.Hallin[as.vector(solve_LSAP(X.DIST)),]
    }
  }
  return(X.rank)
}

discov.Hallin.stat <- function(X, Y, const=1e4, constC=1e8) {
  if (nrow(X) != nrow(Y)) {
    print("Error: sample sizes must agree!")
  } else {
    n = nrow(X)
    X.rank = Hallin.rank(X, const, constC)
    Y.rank = Hallin.rank(Y, const, constC)
    return(n*dcovU(X.rank, Y.rank))
  }
}


############# Nabarun Mrcov ############
######### Computing Rank Energy Statistic #########
computestatisticrdcov=function(x,y,dim1=ncol(x),dim2=ncol(y),n=nrow(x),gridch=halton(n,dim1+dim2),gridch1=as.matrix(gridch[,(1:dim1)]),gridch2=as.matrix(gridch[,((dim1+1):(dim1+dim2))]))
{
  distmat1=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
    distmat1[i,]=apply((x[i,]-t(gridch1)),2,Norm)^2
  assignmentFUN1=solve_LSAP(distmat1)
  assignmentSOL1=cbind(seq_along(assignmentFUN1),assignmentFUN1)
  distmat2=matrix(0,nrow=n,ncol=n)
  for(i in 1:n)
    distmat2[i,]=apply((y[i,]-t(gridch2)),2,Norm)^2
  assignmentFUN2=solve_LSAP(distmat2)
  assignmentSOL2=cbind(seq_along(assignmentFUN2),assignmentFUN2)
  randcovSTAT=dcov.test(gridch1[assignmentSOL1[,2],],gridch2[assignmentSOL2[,2],], R=1)
  return(randcovSTAT$statistic)
}
