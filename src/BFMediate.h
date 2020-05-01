#ifndef __BFMEDIATE_H__
#define __BFMEDIATE_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>

using namespace arma;
using namespace Rcpp;


// shared functions

List runiregGibbs_betafix(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq, double sigmasq, int R, int keep, int nprint, int betafix);

List breg_me(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A);

NumericVector rtrun(NumericVector const& mu, NumericVector const& sigma, NumericVector const& a, NumericVector const& b);

arma::vec dstartoc(arma::vec const& dstar);

double lldstar(arma::vec const& dstar, arma::vec const& y, arma::vec const& mu, double ssq_y_tilde);

double lndMvn(arma::vec const& x, arma::vec const& mu, arma::mat const& rooti);

List dstarRwMetrop(arma::vec const& y, arma::vec const& mu, arma::vec const& olddstar, double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::mat const& rootdi, int ncut, double ssq_y_tilde);

List dstarRwMetrop_1(arma::vec const& y, arma::vec const& mu, arma::vec const& olddstar, double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::mat const& rootdi, int ncut, double ssq_y_tilde);

List breg2(arma::mat const& root, arma::mat const& X, arma::vec const& y, arma::vec const& Abetabar);

double norm_rs(double a, double b);

double half_norm_rs(double a, double b);

double unif_rs(double a, double b);

double exp_rs(double a, double b);

arma::vec rtnm(arma::vec mus, arma::vec sigmas, arma::vec lower, arma::vec upper);

arma::vec rtrunVec(arma::vec const& mu,arma::vec const& sigma, arma::vec const& a, arma::vec const& b);

arma::vec breg1(arma::mat const& root, arma::mat const& X, arma::vec const& y, arma::vec const& Abetabar);

// Main function 1 for Mediation_Ordered_Multi_Merr
List rordprobitGibbs_me(arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad, double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat, int const& Y_ind, int R, int keep, int nprint, arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz);

// Main function 2 for Mediation_Ordered_Multi_Merr
List rordprobitGibbs_me_M(arma::vec const& dep, arma::vec const& beta_2,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad, double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat, int const& Y_ind, int R, int keep, int nprint, arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz);

//FUNCTION TIMING (contained in Timing.cpp)---------------------------------------------------------------
void startMcmcTimer();
void infoMcmcTimer(int rep, int R);
void endMcmcTimer();






#endif