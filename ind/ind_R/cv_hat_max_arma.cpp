#include <RcppArmadillo.h>
#include <math.h>
#include <iomanip>
#include <float.h>
#include <iostream>
#include <algorithm>
#include <random>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;




// [[Rcpp::export]]
arma::mat Max_n(arma::mat X, int p){
  int n = X.n_rows;
  arma::mat A(n, p, arma::fill::zeros);
  arma::uvec row_indices(n, arma::fill::zeros);
  //arma::uvec col_indices = arma::linspace(0,4,5);
  for(int i=0; i<p; i++){
    //A.col(i) = max(X, 1);
    row_indices = index_max(X, 1);
    for(int j=0; j<n; j++){
      A.row(j)(i)=X.row(j)(row_indices(j));
      X.row(j)(row_indices(j))=0;
    } 
  }
  return(A);
  
}


// [[Rcpp::export]]
Rcpp::List cv_arma(arma::mat X, arma::mat Y, arma::mat M, double alpha){
  int n = X.n_rows;
  int p = X.n_cols;
  int q = Y.n_cols;
  int N = M.n_rows;
  int d=min(100, p);
  int tim = p/d;
  int index1; int index2; int index3; int index4;
  double cv_hat_L1; double cv_hat_L3; double cv_hat_L5; 
  arma::mat X_sig;   arma::mat X_block;   arma::mat Y_sig; 
  arma::mat gamma;   arma::mat dat;
  arma::mat dat_new; arma::mat dat_sort;
  arma::mat A;       arma::mat A_sort;
  int l = q/1000;
  if(p <= 500 && q <= 500){
    d = p;  tim = 1; 
  }else if(l > 1){
    d = d/l; tim = p/d;
  };
  A.zeros(N, 5*tim); 
  for(int i=0; i<tim; i++){
    if(i == (tim-1)){
      index1 = i*d;
      index2 = p-1;
    }else{
      index1 = i*d;
      index2 = ((i+1)*d)-1;
    }
    // x rep q 
    X_block = X.cols(index1, index2);
    X_sig = repmat(X_block, q, 1);
    X_sig.reshape(n, d*q);
    // Y rep d
    Y_sig = repmat( Y, 1, d);
    gamma = X_sig % Y_sig;
    dat = gamma - repmat(mean(gamma, 0), n, 1);
    dat_new = abs(M * dat /sqrt(n));
    //dat_sort = sort(dat_new, "descend", 1);
    //1-5 col
    index3 = i*5;
    index4 = ((i+1)*5)-1 ;
    A.cols(index3, index4) = Max_n(dat_new, 5);
  }
  A_sort = Max_n(A, 5);
  arma:: vec A_sotr_1 = A_sort.col(0);
  cv_hat_L1 = Max_n(A_sort.col(0).t(), N*alpha)((N*alpha-1));
  arma::mat C3 =  sum(A_sort.cols(0,2), 1);
  arma::mat C5 =  sum(A_sort.cols(0,4), 1);
  cv_hat_L3 = Max_n(C3.t(), N*alpha)((N*alpha-1));
  cv_hat_L5 = Max_n(C5.t(), N*alpha)((N*alpha-1));
  arma:: vec cv_est = {cv_hat_L1, cv_hat_L3, cv_hat_L5};
  
  return List::create(Named("cv_est") =cv_est);
}






