// #include <Rcpp.h> // 2025-05-27: I'm not an RcppArmadillo coder -GA
#include <iostream>
#include <RcppArmadillo.h>
#include "diagnostics.h"
#include "gibbs_method.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// 2025-05-27: try as best we can to define funcs for taking impact estimates sample by sample

// [[Rcpp::export]]

List impact_cpp_arma(arma::mat betas, arma::vec rhos, arma::sp_mat W) {
  mat direct, indirect, total;
  mat directs, indirects, totals;

  for(int i = 0; i < rhos.n_elem; i++) {
    diagnostic_impacts(betas.row(i), rhos[i], W, direct, indirect, total);

    directs.insert_rows(directs.n_rows, direct);
    indirects.insert_rows(indirects.n_rows, indirect);
    totals.insert_rows(totals.n_rows, total);
  }

  return List::create(Named("direct") = directs,
                      Named("indirect") = indirects,
                      Named("total") = totals);
}

// [[Rcpp::export]]

List impact_Durbin_cpp_arma(arma::mat betas, arma::mat thetas, arma::vec rhos, arma::sp_mat W) {
  mat direct, indirect, total;
  mat directs, indirects, totals;

  for(int i = 0; i < rhos.n_elem; i++) {
    diagnostic_impacts_Durbin(betas.row(i), thetas.row(i), rhos[i], W, direct, indirect, total);

    directs.insert_rows(directs.n_rows, direct);
    indirects.insert_rows(indirects.n_rows, indirect);
    totals.insert_rows(totals.n_rows, total);
  }

  return List::create(Named("direct") = directs,
                      Named("indirect") = indirects,
                      Named("total") = totals);
}
