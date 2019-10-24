#include<RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */


// [[Rcpp::export]]
void armatest(const arma::mat a, arma::colvec vec){
  
  arma::colvec result;
  
  for(int i = 0; i < 5000; ++i){
    arma::colvec result(a*vec);
  }
  
  //https://stackoverflow.com/questions/32579868/return-rcpparmadillo-vector-as-r-vector 
  
  // return result;
  
}

// [[Rcpp::export]]
void armatest2(const arma::sp_mat& a, arma::colvec vec){
  
  // arma::colvec result;
  
  for(int i = 0; i < 5000; ++i){
    arma::colvec result(a*vec);
  }
  
  //https://stackoverflow.com/questions/32579868/return-rcpparmadillo-vector-as-r-vector 
  
  // return result;
  
}


//https://gallery.rcpp.org/articles/armadillo-sparse-matrix-performance/
// [[Rcpp::export]]
arma::sp_mat mult_sp_sp_to_sp(const arma::sp_mat& a, const arma::sp_mat& b) {
  // sparse x sparse -> sparse
  arma::sp_mat result(a * b);
  return result;
}
