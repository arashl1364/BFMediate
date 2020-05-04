#ifndef __BFMEDIATE_H__
#define __BFMEDIATE_H__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
// #include <stdio.h>
// #include <time.h>

using namespace arma;
using namespace Rcpp;

// UtilityFunctions
arma::vec rtrunVec(arma::vec const& mu,arma::vec const& sigma, arma::vec const& a, arma::vec const& b);

arma::vec breg1(arma::mat const& root, arma::mat const& X, arma::vec const& y, arma::vec const& Abetabar);

double lndMvn(arma::vec const& x, arma::vec const& mu, arma::mat const& rooti);

List breg_me(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A);

NumericVector rtrun(NumericVector const& mu, NumericVector const& sigma,
                    NumericVector const& a, NumericVector const& b);

arma::vec dstartoc(arma::vec const& dstar);

double lldstar(arma::vec const& dstar, arma::vec const& y, vec const& mu, double ssq_y_tilde);

List dstarRwMetrop(arma::vec const& y, arma::vec const& mu, arma::vec const& olddstar, double s, arma::mat const& inc_root,
                   arma::vec const& dstarbar, arma::mat const& rootdi, int ncut, double ssq_y_tilde);

List breg2(arma::mat const& root, arma::mat const& X, arma::vec const& y, arma::vec const& Abetabar);

arma::vec rtnm(arma::vec mus, arma::vec sigmas, arma::vec lower, arma::vec upper);

// others
double norm_rs(double a, double b);

double half_norm_rs(double a, double b);

double unif_rs(double a, double b);

double exp_rs(double a, double b);


// Extra functions
List runiregGibbs_betafix(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq,
                          double sigmasq, int R, int keep, int nprint, int betafix);

// Main Functions
// runiregGibbs_rcpp_me.cpp
// [[Rcpp::export]]
List runiregGibbs_rcpp_me(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq,
                          double sigmasq, int R, int keep, int nprint, bool betafix, bool sigmafix, arma::mat betavalue, arma::vec sigmavalue);



// [[Rcpp::export]]
List rordprobitGibbs_me_M_multi_merr_cpp_loop(arma::vec const& dep,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad, arma::mat const& A_2, arma::vec const& betabar_2,
                                              double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                                              int const& Y_ind,
                                              int R, int keep, int nprint);


#endif
