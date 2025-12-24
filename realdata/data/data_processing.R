rm(list=ls()) 
library(tseries)
#install.packages("rstudioapi")  # ????????????
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

name <- c('CS', 'CD', 'CSt', 'Eng', 'Fin', 'HC', 
          'Ind', 'IT', 'Mat', 'RE', 'Uti')
d <- length(name)
de_idx <- NULL



for (k in 1:d) {
  a1 <- paste0("raw_data/68", name[k], ".csv")
  x  <- read.csv(a1)
  ###get stock name
  a <- unique(x$PERMNO)
  num <- rep(0, length(a))
  ###samples of each stocks
  for (i in 1:length(a)) {
    num[i] <- length(which(x$PERMNO==a[i]))
  }
  ind <- which(num!=max(num))
  ####samples are small
  if(length(ind)==0){
    b <- max(num)
    #####4column
    A <- matrix(0, b, length(a)*4)
    X_o <- matrix(0, b, length(a))
    for (j in 1:length(a)) {
      index1 <- (j-1)*4 +1
      index2 <- (j-1)*4 +2
      index3 <- (j-1)*4 +3
      index4 <- j*4
      A[, index1] <- as.character(x$TICKER[which(x$PERMNO==a[j])][1])
      A[, index2] <- x$PERMNO[which(x$PERMNO==a[j])]
      A[, index3] <- x$date[which(x$PERMNO==a[j])]
      A[, index4] <- x$PRC[which(x$PERMNO==a[j])]
      ####stock price
      X_o[, j] <- x$PRC[which(x$PERMNO==a[j])]
    }
  }else{
    ###delete samples are small
    a <- a[-ind] 
    b <- max(num)
    A <- matrix(0, b, length(a)*4)
    X_o <- matrix(0, b, length(a))
    for (j in 1:length(a)) {
      index1 <- (j-1)*4 +1
      index2 <- (j-1)*4 +2
      index3 <- (j-1)*4 +3
      index4 <- j*4
      A[, index1] <- as.character(x$TICKER[which(x$PERMNO==a[j])][1])
      A[, index2] <- x$PERMNO[which(x$PERMNO==a[j])]
      A[, index3] <- x$date[which(x$PERMNO==a[j])]
      A[, index4] <- x$PRC[which(x$PERMNO==a[j])]
      X_o[, j] <- x$PRC[which(x$PERMNO==a[j])]
    }
  }
  a2 <- paste0("raw_data/68new_", name[k], ".csv")
  write.csv(A, a2)
  X <- X_o
  a3 <- paste0("raw_data/68X_", name[k], ".csv")
  write.csv(X, file = a3)
  fun_na <- function(x){
    l <- which(is.na(x)==TRUE)
    if(length(l)>0) {return(length(l))
    }else{return(0)}
  }
  nn <- dim(X)[1]
  x_new <- log(X[2:nn,]/X[1:(nn-1),])
  idx_na <- apply(x_new, 2, fun_na)
  ind_d <- which(idx_na > 0)
  if (length(ind_d) > 0){
    x_new_1 <- x_new[, -ind_d]
    idx_na_1 <- apply(x_new_1, 1, fun_na)
  }else{
    x_new_1 <- x_new
    idx_na_1 <- apply(x_new_1, 1, fun_na)
  }
  #print(which(idx_na_1 > 0))
  ######find all NA-index for all stocks
  if(k == 1){
    de_idx <- which(idx_na_1 > 0)
  }else{
    de_idx <- c(de_idx, which(idx_na_1 > 0))
  }
  #print(k)
  a4 <- paste0("raw_data/68all_", name[k], ".csv")
  write.csv(x_new_1, file = a4)
}
#####find unique NA-index
index_all  <- unique(de_idx)
ind_my <- sort(index_all)

#######final data
for (k in 1: 11) {
  a5 <- paste0("raw_data/68all_", name[k], ".csv")
  x_new_2 <- read.csv(a5)[, -1]
  ###delete all NA with same situation
  if(length(ind_my)>0){
    x_new_3 <- x_new_2[-ind_my, ]
  }else{
    x_new_3 <- x_new_2
  }
  p <- dim(x_new_3)[2]
  res <- NULL
  #####unite root test
  for (i in 1:p) {
    re <- adf.test(x_new_3[,i])
    res[i] <- re$p.value
    
  }
  #####delete non-stationary
  ind_acf <- which(res > 0.01)
  if(length(ind_acf) > 0){
    #print(dim(x_new_3)[2])
    x_new_3 <- x_new_3[, -ind_acf]
    #print(dim(x_new_3)[2])
  }
  a6 <- paste0("processed_data/68return_", name[k], ".csv")
  write.csv(x_new_3, file = a6)
}


########## second pierd2020-2022

for (k in 1:d) {
  a1 <- paste0("./raw_data/02", name[k], ".csv")
  x  <- read.csv(a1)
  ###get stock name
  a <- unique(x$PERMNO)
  num <- rep(0, length(a))
  ###samples of each stocks
  for (i in 1:length(a)) {
    num[i] <- length(which(x$PERMNO==a[i]))
  }
  ind <- which(num!=max(num))
  ####samples are small
  if(length(ind)==0){
    b <- max(num)
    #####4column
    A <- matrix(0, b, length(a)*4)
    X_o <- matrix(0, b, length(a))
    for (j in 1:length(a)) {
      index1 <- (j-1)*4 +1
      index2 <- (j-1)*4 +2
      index3 <- (j-1)*4 +3
      index4 <- j*4
      A[, index1] <- as.character(x$TICKER[which(x$PERMNO==a[j])][1])
      A[, index2] <- x$PERMNO[which(x$PERMNO==a[j])]
      A[, index3] <- x$date[which(x$PERMNO==a[j])]
      A[, index4] <- x$PRC[which(x$PERMNO==a[j])]
      ####stock price
      X_o[, j] <- x$PRC[which(x$PERMNO==a[j])]
    }
  }else{
    ###delete samples are small
    a <- a[-ind] 
    b <- max(num)
    A <- matrix(0, b, length(a)*4)
    X_o <- matrix(0, b, length(a))
    for (j in 1:length(a)) {
      index1 <- (j-1)*4 +1
      index2 <- (j-1)*4 +2
      index3 <- (j-1)*4 +3
      index4 <- j*4
      A[, index1] <- as.character(x$TICKER[which(x$PERMNO==a[j])][1])
      A[, index2] <- x$PERMNO[which(x$PERMNO==a[j])]
      A[, index3] <- x$date[which(x$PERMNO==a[j])]
      A[, index4] <- x$PRC[which(x$PERMNO==a[j])]
      X_o[, j] <- x$PRC[which(x$PERMNO==a[j])]
    }
  }
  a2 <- paste0("raw_data/02new_", name[k], ".csv")
  write.csv(A, a2)
  X <- X_o
  a3 <- paste0("raw_data/02X_", name[k], ".csv")
  write.csv(X, file = a3)
  fun_na <- function(x){
    l <- which(is.na(x)==TRUE)
    if(length(l)>0) {return(length(l))
    }else{return(0)}
  }
  nn <- dim(X)[1]
  x_new <- log(X[2:nn,]/X[1:(nn-1),])
  idx_na <- apply(x_new, 2, fun_na)
  ind_d <- which(idx_na > 0)
  if (length(ind_d) > 0){
    x_new_1 <- x_new[, -ind_d]
    idx_na_1 <- apply(x_new_1, 1, fun_na)
  }else{
    x_new_1 <- x_new
    idx_na_1 <- apply(x_new_1, 1, fun_na)
  }
  #print(which(idx_na_1 > 0))
  ######find all NA-index for all stocks
  if(k == 1){
    de_idx <- which(idx_na_1 > 0)
  }else{
    de_idx <- c(de_idx, which(idx_na_1 > 0))
  }
  #print(k)
  a4 <- paste0("raw_data/02all_", name[k], ".csv")
  write.csv(x_new_1, file = a4)
}
#####find unique NA-index
index_all  <- unique(de_idx)
ind_my <- sort(index_all)

#######final data
for (k in 1: 11) {
  a5 <- paste0("raw_data/02all_", name[k], ".csv")
  x_new_2 <- read.csv(a5)[, -1]
  ###delete all NA with same situation
  if(length(ind_my)>0){
    x_new_3 <- x_new_2[-ind_my, ]
  }else{
    x_new_3 <- x_new_2
  }
  p <- dim(x_new_3)[2]
  res <- NULL
  #####unite root test
  for (i in 1:p) {
    re <- adf.test(x_new_3[,i])
    res[i] <- re$p.value
    
  }
  #####delete non-stationary
  ind_acf <- which(res > 0.01)
  if(length(ind_acf) > 0){
    #print(dim(x_new_3)[2])
    x_new_3 <- x_new_3[, -ind_acf]
    #print(dim(x_new_3)[2])
  }
  a6 <- paste0("processed_data/02return_", name[k], ".csv")
  write.csv(x_new_3, file = a6)
}
