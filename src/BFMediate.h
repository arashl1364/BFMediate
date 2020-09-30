#ifndef __BFMEDIATE_H__
#define __BFMEDIATE_H__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <time.h>

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

List dstarRwMetrop2(arma::vec const& y, arma::vec const& mu, arma::vec const& olddstar, double s, arma::mat const& inc_root,
                    arma::vec const& dstarbar, double oldll, arma::mat const& rootdi, int ncut, double ssq_y_tilde);

List dstarRwMetrop3(arma::vec const& y, arma::vec const& mu, arma::vec const& olddstar, double s, arma::mat const& inc_root,
                    arma::vec const& dstarbar, arma::mat const& rootdi, int ncut, double ssq_y_tilde);

List breg2(arma::mat const& root, arma::mat const& X, arma::vec const& y, arma::vec const& Abetabar);

arma::vec rtnm(arma::vec mus, arma::vec sigmas, arma::vec lower, arma::vec upper);

List breg2(arma::mat const& root, arma::mat const& X, arma::vec const& y, arma::vec const& Abetabar);

List runiregG(arma::vec const& y, arma::mat const& X, arma::mat const& XpX, arma::vec const& Xpy, double sigmasq, arma::mat const& A,
              arma::vec const& Abetabar, double nu, double ssq);
// others
double norm_rs(double a, double b);

double half_norm_rs(double a, double b);

double unif_rs(double a, double b);

double exp_rs(double a, double b);


// Truncated normal draw
double dexpr(double const& a);

double invCdfNorm(double const& a);

double dnr(double const& a);

double trunNormBelow(double const& a);

double trunNorm(double mu,double sig, double trunpt, int above);

arma::vec trunNorm_vec(arma::vec const& mu, arma::vec const& sig, arma::vec const& trunpt, arma::vec const& above);



// Extra functions
List runiregGibbs_betafix(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq,
                          double sigmasq, int R, int keep, int nprint, int betafix);

// List rordprobitGibbs_me
List YSampler(arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
                        double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                        int const& Y_ind,
                        int R, int keep, int nprint,
                        arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz);

List MSampler(arma::vec const& dep, arma::vec const& beta_2,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
                          double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                          int const& Y_ind,
                          int R, int keep, int nprint,
                          arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz);

// List rordprobitGibbs_me_M(arma::vec const& dep, arma::vec const& beta_2,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
//                           double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
//                           int const& Y_ind,
//                           int R, int keep, int nprint,
//                           arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz);


//FUNCTION TIMING (contained in TimerFunctions.cpp)---------------------------------------------------------------
void startMcmcTimer();
void infoMcmcTimer(int rep, int R);
void endMcmcTimer();

// Main Functions
// runiregGibbs_rcpp_me.cpp
// [[Rcpp::export]]
List runiregGibbs_rcpp_me(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq,
                          double sigmasq, int R, int keep, int nprint, bool betafix, bool sigmafix, arma::mat betavalue, arma::vec sigmavalue);



// //  MeasurementMCatUnitCpp.cpp
// // [[Rcpp::export]]
// List MeasurementMCatUnitCpp(arma::vec const& dep,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad, arma::mat const& A_2, arma::vec const& betabar_2,
//                                               double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
//                                               int const& Y_ind,
//                                               int R, int keep, int nprint);

//  MeasurementMCatCpp.cpp
// [[Rcpp::export]]
List MeasurementMCatCpp(arma::vec const& dep,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad, arma::mat const& A_2, arma::vec const& betabar_2,
                            double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                            int const& Y_ind,
                            int R, int keep, int nprint);

// MeasurementYCatCpp.cpp
// [[Rcpp::export]]
List MeasurementYCatCpp(arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
                                            double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                                            int const& Y_ind,
                                            int R, int keep, int nprint);
// MeasurementMYCatCpp.cpp
// [[Rcpp::export]]
List MeasurementMYCatCpp(arma::mat const& X, arma::mat const& m_star, arma::mat const& y_star, int k_M, int k_Y, int M_ind, int Y_ind,     //data
                                      arma::mat const& A_M, arma::vec const& betabar, arma::mat const& Ad_M, double s_M, arma::mat const& inc_root_M, arma::vec const& dstarbar_M, arma::vec const& betahat,               //priors_M
                                      arma::mat const& A_Y, arma::vec const& beta_2_bar, arma::mat const& Ad_Y, double s_Y, arma::mat const& inc_root_Y, arma::vec const& dstarbar_Y, arma::vec const& beta_2_hat,         //priors_Y
                                      int R, int keep, int nprint);
// RuniregGibbsMultiCpp.cpp
// [[Rcpp::export]]
List RuniregGibbsMultiCpp(arma::mat const& y, arma::mat const& M, arma::vec const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq,
                             arma::vec sigmasq, int R, int keep, int nprint, bool betafix, bool sigmafix, arma::mat betavalue, arma::vec sigmavalue);



#endif
