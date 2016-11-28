// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-

#include <memory>
#include "RcppArmadillo.h"

#include "kmpp.hpp"

#include "mixture.hpp"
#include "reducedDimARC.hpp"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
Rcpp::List jkmeansEM(const arma::mat& y, int k, const arma::mat& meansIni,
                     int steps = 1000, double tol = 1E-8, bool fixW = false,
                     bool useKmeansIni = false, double sigma2_ini = 0.1) {
  Mixture mix(y, k);

  mix.initialize(meansIni, useKmeansIni, fixW, sigma2_ini);

  mix.runEM(steps, tol);

  mix.sortBy1stDinMu();

  return Rcpp::List::create(
      Rcpp::Named("mu") = mix.mu, Rcpp::Named("w") = mix.w,
      Rcpp::Named("Sigma") = mix.Sigma, Rcpp::Named("zeta") = mix.zeta,
      Rcpp::Named("M") = mix.clusteringMAP());
}

// [[Rcpp::export]]
Rcpp::List rDARC(const arma::mat& Y, int d, int k, const arma::mat& meansIni,
                 int steps = 1000, double tol = 1E-8, bool randomStart = true,
                 bool fixW = false, double sigma2_ini = 0.1, int ver = 1) {
  reducedDimARC rdARC(Y, d, k, randomStart, sigma2_ini, ver, fixW);
  rdARC.runEM(steps);

  return Rcpp::List::create(
      Rcpp::Named("mu") = rdARC.mix->mu, Rcpp::Named("w") = rdARC.mix->w,
      Rcpp::Named("Sigma") = rdARC.mix->Sigma,
      Rcpp::Named("zeta") = rdARC.mix->zeta, Rcpp::Named("X") = rdARC.X,
      Rcpp::Named("EV") = rdARC.EV,
      Rcpp::Named("M") = rdARC.mix->clusteringMAP()

          );
}

// [[Rcpp::export]]
Rcpp::List kmpp(const arma::mat& Y, int k) {
  return Rcpp::List::create(Rcpp::Named("mu") = kmeanspp(Y, k));
}
// // [[Rcpp::export]]
// Rcpp::List jkmeansEMBatch(const arma::cube& y, int k, int j, int steps =
// 1000,
//                           double tol = 1E-8, bool fixW = true,
//                           bool flexJ = false, double zetaTrunc = 0.01,
//                           bool useKmeansIni = true,
//                           const arma::mat& meansIni = 0,
//                           double sigma2_ini = 0.1, bool normalizeZeta =
//                           false) {
//   if (j > k) {
//     throw std::range_error("j needs be no bigger than k");
//   }

//   int n = y.n_rows;
//   int p = y.n_cols;
//   int batchN = y.n_slices;

//   cube mu(k, p, batchN);
//   mat w(k, batchN);
//   vec sigma2(batchN);
//   cube zeta(n, k, batchN);
//   umat M(n, batchN);

// #pragma omp parallel for
//   for (int i = 0; i < batchN; ++i) {
//     mat localY = y.slice(i);
//     Mixture mix(localY, k, j);
//     mix.initialize(meansIni, useKmeansIni, fixW, flexJ, zetaTrunc,
//     sigma2_ini,
//                    normalizeZeta);
//     mix.runEM(steps, tol);

//     mix.sortBy1stDinMu();

//     mu.slice(i) = mix.mu;
//     w.col(i) = mix.w;
//     sigma2(i) = mix.sigma2;
//     zeta.slice(i) = mix.zeta;
//     M.col(i) = mix.clusteringMAP();
//   }

//   return Rcpp::List::create(Rcpp::Named("mu") = mu, Rcpp::Named("w") = w,
//                             Rcpp::Named("sigma2") = sigma2,
//                             Rcpp::Named("zeta") = zeta, Rcpp::Named("M") =
//                             M);
// }