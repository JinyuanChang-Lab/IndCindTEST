###packages
library(clue)
library(adagio)
library(energy)
library(pracma)
library(randtoolbox)
library(RcppArmadillo)



#---------------------------------------------------------- The competing methods for independence test --------------------------

######################     Zhu et al. (2017) -Pcor- provided by authors   ####################
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



################ Chakraborty and Zhang (2019) -Jdcov_R- provided in the supplementary material of Chakraborty and Zhang (2019).  ##########
#####https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1513364#supplemental-material-section
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



###################  Shi et al. (2022) -Hallin- provided in tha supplementary material of Shi et al. (2022) ######
#####https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1782223#supplemental-material-section
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


############# Deb and Sen (2023) -mrdCov- provided in the supplementary material of Deb and Sen (2023) #####
############ https://github.com/NabarunD
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




#--------------------------------- Implement independence test of x and y using the proposed method and competing methods -------

Gauss_T <- function(n, p, q, dep, alpha, rep.sim){
  
  sA <- rep(0, 12)                                                           # generate a 12-by-1 vector  
  sB <- rep(0, 12)                                                           # generate a 12-by-1 vector
  Jd <- rep(0, 3)                                                            # generate a 3-by-1 vector
  
  BstrpTimes <- 500                                                          # Number of bootstrap/permutation iterations
  Bstrp <-Bstrp_gaussian <-  matrix(0,BstrpTimes,1)                          # generate a BstrpTimes-by-1 matrix
  Bs_han <- Bs_han_gaussian <- matrix(0,BstrpTimes,1)                        # generate a BstrpTimes-by-1 matrix
  Bs_Mr <- Bs_Mr_gaussian <- matrix(0,BstrpTimes,1)                          # generate a BstrpTimes-by-1 matrix
  
  set.seed(rep.sim+1000)
  dat <- Data1(n, p, q, dep)                                                 # Generate random data by the function Data1
  x <- dat$x                                                                 # X matrix (n-by-p)
  y <- dat$y                                                                 # Y matrix (n-by-q)
  
  
 
  
  #coordinatewise Gasuusianization
  x_norm <- data_transform(x)$x_trs_my                                        
  y_norm <- data_transform(y)$x_trs_my                                        
  
  ############# the proposed independence test ############
  seed_my = rep.sim+100
  res_my <- Ind_Gtest(x, y, alpha, seed_my,  "all")
  sA[1] <- res_my$Gaussian$reject      # Gaussian multiplier
  sA[2] <- res_my$Mammen$reject        # Mammen's multiplier 
  sA[3] <- res_my$Rademacher$reject    # Rademacher multiplier
  
  
  #################### Pcor ########################
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
  
  ######## Pcor--original data
  set.seed(rep.sim+1000)
  pcov_inv <- try(pcov_fun(x,y,BstrpTimes),  silent=TRUE) 
  if ('try-error' %in% class(pcov_inv)) {
    pVal <- 8888      # indicates failure; downstream code must check and skip  
  }else{ 
    pVal = pcov_inv
  }
  
  ######## Pcor--coordinatewise Gasuusianization 
  set.seed(rep.sim+1000)
  pcov_inv_gaussian <- try(pcov_fun(x_norm, y_norm, BstrpTimes),  silent=TRUE) 
  if ('try-error' %in% class(pcov_inv_gaussian)) {
    pVal_gaussian <- 8888  # indicates failure; downstream code must check and skip
  }else{ 
    pVal_gaussian = pcov_inv_gaussian
  }
  
  
  ############################### rdCov ########################
  dr_fun <- function(x,y, BstrpTimes){
    Dx = as.matrix(dist((x)), diag = TRUE, upper = TRUE)
    Dy = as.matrix(dist((y)), diag = TRUE, upper = TRUE)
    hhg = hhg.test(Dx, Dy, nr.perm =BstrpTimes)
    drpVal =hhg$perm.pval.hhg.sc;
    return(drpVal)
  } 
  
  ######## rdCov--original data
  set.seed(rep.sim+1000)
  dr_inv <- try(dr_fun(x,y, BstrpTimes),  silent=TRUE) 
  if ('try-error' %in% class(dr_inv)) {
    drpVal <- 8888   # indicates failure; downstream code must check and skip
  }else{ 
    drpVal <- ifelse(is.nan(dr_inv), 9999, dr_inv)
  }
  
  ######## rdCov--coordinatewise Gasuusianization 
  set.seed(rep.sim+1000)
  dr_inv_gaussian <- try(dr_fun(x_norm,y_norm, BstrpTimes),  silent=TRUE) 
  if ('try-error' %in% class(dr_inv_gaussian)) {
    drpVal_gaussian <- 8888    # indicates failure; downstream code must check and skip
  }else{ 
    drpVal_gaussian <- ifelse(is.nan(dr_inv_gaussian), 9999, dr_inv_gaussian)
  }
  
  
  ######################## dCor ########################
  ###### dCor--original data
  set.seed(rep.sim+1000)
  dcor_inv <- try(dcorT.test(x,y)$p.value,  silent=TRUE)
  if ('try-error' %in% class(dcor_inv)) {
    dcortpVal <- 8888    # indicates failure; downstream code must check and skip
  }else{ 
    dcortpVal <- ifelse(is.nan(dcor_inv), 9999, dcor_inv)
  }
  
  ###### dCor--coordinatewise Gasuusianization
  set.seed(rep.sim+1000)
  dcor_inv_gaussian <- try(dcorT.test(x_norm,y_norm)$p.value,  silent=TRUE)
  if ('try-error' %in% class(dcor_inv_gaussian)) {
    dcortpVal_gaussian <- 8888      # indicates failure; downstream code must check and skip
  }else{ 
    dcortpVal_gaussian <- ifelse(is.nan(dcor_inv_gaussian), 9999, dcor_inv_gaussian)
  }
  
  ############### dHSIC #####################
  
  ########### dHSIC--original data
  set.seed(rep.sim+1000)
  dhsic_inv <- try(dhsic.test(x, y)$p.value,  silent=TRUE)
  if ('try-error' %in% class(dhsic_inv)) {
    dhsicpVal <- 8888    # indicates failure; downstream code must check and skip
  }else{ 
    dhsicpVal <- ifelse(is.nan(dhsic_inv)|is.na(dhsic_inv), 9999, dhsic_inv)
  }
  
  ########### dHSIC--coordinatewise Gasuusianization
  set.seed(rep.sim+1000)
  dhsic_inv_gaussian <- try(dhsic.test(x_norm, y_norm)$p.value,  silent=TRUE)
  if ('try-error' %in% class(dhsic_inv_gaussian)) {
    dhsicpVal_gaussian <- 8888   # indicates failure; downstream code must check and skip
  }else{ 
    dhsicpVal_gaussian <- ifelse(is.nan(dhsic_inv_gaussian)|is.na(dhsic_inv_gaussian), 9999, dhsic_inv_gaussian)
  }
  
  
  
  z <- cbind(x,y)
  d <- dim(z)[2]
  z_norm <- cbind(x_norm,y_norm)
  
  ##################### JdCov_R ########################
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
    stat.p <- stat.pr <- stat.ps <- rep(0, Mp)
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
  
  ###### JdCov_R--original data  
  set.seed(rep.sim+1000)
  jdcov_inv <- try(sta13(z, Mp, n, d, alpha),  silent=TRUE) 
  if ('try-error' %in% class(jdcov_inv)) {
    Jd <- rep(8888, 3)  # indicates failure; downstream code must check and skip
  } else{ 
    Jd <- jdcov_inv
  }
  
  
  ########### JdCov_R--coordinatewise Gasuusianization
  set.seed(rep.sim+1000)
  jdcov_inv_gaussian <- try(sta13(z_norm, Mp, n, d, alpha),  silent=TRUE) 
  if ('try-error' %in% class(jdcov_inv_gaussian)) {
    Jd_gaussian <- rep(8888, 3) # indicates failure; downstream code must check and skip
  } else{ 
    Jd_gaussian <- jdcov_inv_gaussian
  }
  
  
  ################# GdCov #############################
  
  ################## GdCov--original data  
  set.seed(rep.sim+1000)
  m.inv <- try(mdm_test(z, c(p,q))$pval, silent=TRUE)  
  if ('try-error' %in% class(m.inv)) {
    dmdmpVal <- 8888   # indicates failure; downstream code must check and skip
  } else{ 
    dmdmpVal <- ifelse(is.nan(m.inv)|is.na(m.inv), 9999, m.inv)
  }
  
  ########### GdCov--coordinatewise Gasuusianization
  set.seed(rep.sim+1000)
  m.inv_gaussian <- try(mdm_test(z_norm, c(p,q))$pval, silent=TRUE)  
  if ('try-error' %in% class(m.inv_gaussian)) {
    dmdmpVal_gaussian <- 8888  # indicates failure; downstream code must check and skip
  } else{ 
    dmdmpVal_gaussian <- ifelse(is.nan(m.inv_gaussian)|is.na(m.inv_gaussian), 9999, m.inv_gaussian)
  }
  
  #################### Hallin ##############################
  
  ######### Hallin--original data  
  set.seed(rep.sim+1000)
  had.inv <- try({ 
    H_dcovt <- discov.Hallin.stat(x,y)
    for (j in 1:BstrpTimes){
      O1=order(runif(n,0,1))
      YBstrp1 = y[O1,]
      Bs_han[j]  = discov.Hallin.stat(x,YBstrp1)
    } 
  } , silent=TRUE)
  
  if ('try-error' %in% class(had.inv)) {
    HaVal <- 8888  # indicates failure; downstream code must check and skip
  }else{ 
    HaVal <- mean(Bs_han >= H_dcovt)
  }
  
  ######### Hallin--coordinatewise Gasuusianization
  set.seed(rep.sim+1000)
  had.inv_gaussian <- try({ 
    H_dcovt_gaussian <- discov.Hallin.stat(x_norm,y_norm)
    for (j in 1:BstrpTimes){
      O1=order(runif(n,0,1))
      YBstrp1_gaussian = y_norm[O1,]
      Bs_han_gaussian[j]  = discov.Hallin.stat(x_norm,YBstrp1_gaussian)
    } 
  } , silent=TRUE)
  
  if ('try-error' %in% class(had.inv_gaussian)) {
    HaVal_gaussian <- 8888  # indicates failure; downstream code must check and skip
  }else{ 
    HaVal_gaussian <- mean(Bs_han_gaussian >= H_dcovt_gaussian)
  }
  
  
  ####################### mrdCov ###################################
  
  ###### mrdCov--original data
  set.seed(rep.sim+1000)
  mrn.inv <- try( { 
    Mrcov <- computestatisticrdcov(x, y)
    for (j in 1:BstrpTimes){
      O2=order(runif(n,0,1));
      YBstrp2 = y[O2,];
      Bs_Mr[j]  = computestatisticrdcov(x, YBstrp2);
    }
  } , silent=TRUE)
  
  if ('try-error' %in% class(mrn.inv)) {
    MrVal <- 8888  # indicates failure; downstream code must check and skip
  }else{ 
    MrVal <- mean(Bs_Mr >= Mrcov)
  }
  
  ###### mrdCov--coordinatewise Gasuusianization
  set.seed(rep.sim+1000)
  mrn.inv_gaussian <- try( { 
    Mrcov <- computestatisticrdcov(x_norm, y_norm)
    for (j in 1:BstrpTimes){
      O2=order(runif(n,0,1));
      YBstrp2_gaussian = y_norm[O2,];
      Bs_Mr_gaussian[j]  = computestatisticrdcov(x_norm, YBstrp2_gaussian);
    }
  } , silent=TRUE)
  
  if ('try-error' %in% class(mrn.inv_gaussian)) {
    MrVal_gaussian <- 8888  # indicates failure; downstream code must check and skip
  }else{ 
    MrVal_gaussian <- mean(Bs_Mr >= Mrcov)
  }
  
  resn <- c(pVal, drpVal, dcortpVal, dhsicpVal,
            Jd, dmdmpVal, HaVal, MrVal, sA, sB)
  resn_gaussian <- c(pVal_gaussian, drpVal_gaussian, dcortpVal_gaussian, dhsicpVal_gaussian,
                     Jd_gaussian, dmdmpVal_gaussian, HaVal_gaussian, MrVal_gaussian, sA, sB)
  result <- list(resn = resn, resn_gaussian = resn_gaussian)
  return(result)
}
