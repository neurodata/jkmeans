// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-

#include "RcppArmadillo.h"
#include "mixture.hpp"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
Rcpp::List jkmeans(const arma::mat& y, int k, int j, int steps = 1000) {
  Mixture mix(y, k, j);

    if (j > k) {          
        throw std::range_error("j needs be no bigger than k");
    }

  mix.run(steps);

  return Rcpp::List::create(Rcpp::Named("mu") = mix.mu,
                            Rcpp::Named("w") = mix.w,
                            Rcpp::Named("sigma2") = mix.sigma2);
}
