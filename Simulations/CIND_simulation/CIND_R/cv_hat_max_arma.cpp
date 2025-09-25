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
arma::mat Max_n(arma::mat X, int p){             // Compute the top-p values in each row of the matrix X
  int n = X.n_rows;                              // number of X rows
  arma::mat A(n, p, arma::fill::zeros);          // Generate an n-by-p matrix
  arma::uvec row_indices(n, arma::fill::zeros);  // Generate an n-by-1 vector
  
  
  for(int i=0; i<p; i++){
    //A.col(i) = max(X, 1);
    row_indices = index_max(X, 1);               // Return the index of the maximum value for each row of X
    for(int j=0; j<n; j++){
      A.row(j)(i)=X.row(j)(row_indices(j));      // A[j,i] = the i-th largest entry of X's j-th row
      X.row(j)(row_indices(j))=0;                // set the i-th largest entry of X's j-th row to zero
    }                                 
  }
  return(A);                                    
  
}




// [[Rcpp::export]]
Rcpp::List cv_arma(arma::mat X, arma::mat Y, arma::mat M, double alpha){   // Compute critical values via multiplier bootstrap; M: multiplier matrix (N-by-n)
  int n = X.n_rows;                                                        // sample size
  int p = X.n_cols;                                                        // number of X columns 
  int q = Y.n_cols;                                                        // number of Y columns 
  int N = M.n_rows;                                                        // number of bootstrap replications using the multiplier bootstrap method                     
  int index1; int index2; int index3; int index4;                           
  double cv_hat_L1; double cv_hat_L3; double cv_hat_L5;                    
  arma::mat X_sig;   arma::mat X_block;   arma::mat Y_sig;                  
  arma::mat gamma;   arma::mat dat;                                         
  arma::mat dat_new; arma::mat dat_sort;                                   
  arma::mat A;       arma::mat A_sort;                                      
  int d=min(100, p);                                                       // We need to split X into column-wise blocks to reduce memory usage when p and q are very large. d is the block size (minimum 100)
  int tim = p/d;                                                           // number of blocks
  int l = q/1000;                                                          // scale factor to further shrink block size if q is very large
  if(p <= 500 && q <= 500){                                                // If p and q are <= 500, no block partitioning is applied.
    d = p;  tim = 1; 
  }else if(l > 1){                                                         // If q > 1000, reduce block size accordingly.
    d = d/l; tim = p/d;                                                    
  };
  A.zeros(N, 5*tim);                                                       // Generate an n-by-(5*tim) matrix 
  for(int i=0; i<tim; i++){                                                // loop over X column-wise blocks
    // Determine the column range of the i-th block of X; index1 is the starting column, index2 is the ending column
    if(i == (tim-1)){                                                      
      index1 = i*d;                                                        // the last block
      index2 = p-1;                                                        
    }else{
      index1 = i*d;                                                        // the first (tim-1) blocks
      index2 = ((i+1)*d)-1;                                                
    }
    // 
    X_block = X.cols(index1, index2);                                      // extract the i-th bloch of X (n-by-d matrix)
    X_sig = repmat(X_block, q, 1);                                         // replicate X_block q times ((n*q)-by-d matrix)
    X_sig.reshape(n, d*q);                                                 // reshape X_sig to an n-by-(d*q) matrix
    // 
    Y_sig = repmat( Y, 1, d);                                              // replicate Y d times (n-by-(d*q) matrix)
    // 
    gamma = X_sig % Y_sig;                                                 // Hadamard products for X_sig and Y_sig
    dat = gamma - repmat(mean(gamma, 0), n, 1);                            // Subtract the column mean from each column of gamma
    dat_new = abs(M * dat /sqrt(n));                                       // multiplier transform (N-by-(d*q) matrix)
    //dat_sort = sort(dat_new, "descend", 1);                               
    
    index3 = i*5;                                                           
    index4 = ((i+1)*5)-1 ;                                                  
    A.cols(index3, index4) = Max_n(dat_new, 5);                            // A[index3, index4] stores the top-5 values in each row of dat_new 
  }
  A_sort = Max_n(A, 5);                                                    // take the top-5 values in each row of A
  // L1: use the largest column (top-1) of A_sort
  arma:: vec A_sotr_1 = A_sort.col(0);                                     // rowwise largest value of A_sort (N-by-1 vector) 
  cv_hat_L1 = Max_n(A_sort.col(0).t(), N*alpha)((N*alpha-1));              // [N*alpha]-th largest value of A_sotr_1 
  // L3/L5: rowwise sums of top-3 / top-5 
  arma::mat C3 =  sum(A_sort.cols(0,2), 1);                                // calculate the sum of the top-3 values in each row of A_sort
  arma::mat C5 =  sum(A_sort.cols(0,4), 1);                                // calculate the sum of the top-5 values in each row of A_sort
  cv_hat_L3 = Max_n(C3.t(), N*alpha)((N*alpha-1));                         // [N*alpha]-th largest value of C3
  cv_hat_L5 = Max_n(C5.t(), N*alpha)((N*alpha-1));                         // [N*alpha]-th largest value of C5
  
  // calculate p value
  double ts = arma::abs((X.t() * Y) / std::sqrt(static_cast<double>(n))).max();                  // calculate the test statistic
  double prop = arma::accu(A_sotr_1 > ts) / static_cast<double>(N);                              // calculate the proportion of elements in A_sotr_1 greater than ts   
  
  
  return Rcpp::List::create(  // return as an R list
    Rcpp::Named("ts")      = ts,
    Rcpp::Named("p_value") = prop,
    Rcpp::Named("cv_e")    = cv_hat_L1 
  );                          
}







