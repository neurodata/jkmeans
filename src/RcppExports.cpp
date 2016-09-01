// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// jkmeansEM
Rcpp::List jkmeansEM(const arma::mat& y, int k, int j, int steps, double tol, bool fixW, bool useKmeansIni, const arma::mat& meansIni);
RcppExport SEXP jkmeans_jkmeansEM(SEXP ySEXP, SEXP kSEXP, SEXP jSEXP, SEXP stepsSEXP, SEXP tolSEXP, SEXP fixWSEXP, SEXP useKmeansIniSEXP, SEXP meansIniSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type fixW(fixWSEXP);
    Rcpp::traits::input_parameter< bool >::type useKmeansIni(useKmeansIniSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type meansIni(meansIniSEXP);
    __result = Rcpp::wrap(jkmeansEM(y, k, j, steps, tol, fixW, useKmeansIni, meansIni));
    return __result;
END_RCPP
}
// jkmeansEMBatch
Rcpp::List jkmeansEMBatch(const arma::cube& y, int k, int j, int steps, double tol, bool fixW, bool useKmeansIni, const arma::mat& meansIni);
RcppExport SEXP jkmeans_jkmeansEMBatch(SEXP ySEXP, SEXP kSEXP, SEXP jSEXP, SEXP stepsSEXP, SEXP tolSEXP, SEXP fixWSEXP, SEXP useKmeansIniSEXP, SEXP meansIniSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::cube& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type fixW(fixWSEXP);
    Rcpp::traits::input_parameter< bool >::type useKmeansIni(useKmeansIniSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type meansIni(meansIniSEXP);
    __result = Rcpp::wrap(jkmeansEMBatch(y, k, j, steps, tol, fixW, useKmeansIni, meansIni));
    return __result;
END_RCPP
}
