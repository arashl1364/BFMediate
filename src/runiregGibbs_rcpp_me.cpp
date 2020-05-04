#include "BFMediate.h"
// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// using namespace arma; // use the Armadillo library for matrix computations
// using namespace Rcpp;



List runiregGibbs_rcpp_me(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq,
                               double sigmasq, int R, int keep, int nprint, bool betafix, bool sigmafix, arma::mat betavalue, arma::vec sigmavalue) {

  // Keunwoo Kim 09/09/2014
  //edited by Arash Laghaie 13/02/2018 for estimating the mediation model with measurement error

  // Purpose: perform iid draws from posterior of regression model using conjugate prior

  // Arguments:
  //  y,X
  //  betabar,A      prior mean, prior precision
  //  nu, ssq        prior on sigmasq
  //  R number of draws
  //  keep thinning parameter

  // Output: list of beta, sigmasq

  // Model:
  //  y = Xbeta + e  e ~N(0,sigmasq)
  //  y is n x 1
  //  X is n x k
  //  beta is k x 1 vector of coefficients

  // Prior:
  //  beta ~ N(betabar,sigmasq*A^-1)
  //  sigmasq ~ (nu*ssq)/chisq_nu
  //
  int mkeep;
  double s;
  arma::mat RA, W, IR;
  arma::vec z, btilde;
  //vec beta(2);
  arma::vec beta(X.n_cols);

  int nvar = X.n_cols;
  int nobs = y.size();

  arma::vec sigmasqdraw(R/keep);
  arma::mat betadraw(R/keep, nvar);
  arma::cube IRdraw(nvar,nvar,R/keep);
  arma::mat btildedraw(R/keep,nvar);
  arma::vec nudraw(R/keep);  //nu1 posterior degrees of freedom
  arma::vec Sdraw(R/keep); //S1 posterior scale

  arma::mat XpX = trans(X)*X;
  arma::vec Xpy = trans(X)*y;

  arma::vec Abetabar = A*betabar;


  // if (nprint>0) startMcmcTimer();

  // beta(0)=0;    //we know apriori that the  intercept is equal to 0 and coefficient of M is equal to 1
  // beta(1)=1;


  for (int rep=0; rep<R; rep++){

    if(betafix == true) beta = trans(betavalue.row(rep));
    if(sigmafix == true) sigmasq = sigmavalue(rep);

    //first draw beta | sigmasq
    IR = solve(trimatu(chol(XpX/sigmasq+A)), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
    btilde = (IR*trans(IR)) * (Xpy/sigmasq+Abetabar);
    IRdraw.slice(rep) = IR*trans(IR); //conv_to<uword>::from(IR);
    btildedraw.row(rep) = trans(btilde);

    if(betafix == false) beta = btilde + IR*vec(rnorm(nvar));

    //now draw sigmasq | beta
    s = sum(square(y-X*beta));
    nudraw(rep) = nu+nobs;
    Sdraw(rep) = (nu*ssq+s)/(nu+nobs);

    if(sigmafix == false) sigmasq = (nu*ssq+s) / rchisq(1,nu+nobs)[0]; //rchisq returns a vectorized object, so using [0] allows for the conversion to double

    //print time to completion and draw # every nprint'th draw
    // if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);

    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw(mkeep-1, span::all) = trans(beta);
      sigmasqdraw[mkeep-1] = sigmasq;
    }
  }

  // if (nprint>0) endMcmcTimer();

  return List::create(
    Named("betadraw") = betadraw,
    Named("sigmasqdraw") = NumericVector(sigmasqdraw.begin(),sigmasqdraw.end()),
    //moments for computing Bayes factor
    Named("mubeta") = btildedraw,
    Named("IR")  = IRdraw,
    Named("nu") = nudraw,
    Named("S") = Sdraw);

}



// //The functions below are used to print the output from MCMC draws for many of the bayesm functions
// time_t itime;
// char buf[100];
//
// //[[Rcpp::export]]
// void startMcmcTimer() {
//   itime = time(NULL);
//   Rcout << " MCMC Iteration (est time to end - min) \n";
// }
//
// //[[Rcpp::export]]
// void infoMcmcTimer(int rep, int R) {
//   time_t ctime = time(NULL);
//   char buf[32];
//
//   double timetoend = difftime(ctime, itime) / 60.0 * (R - rep - 1) / (rep+1);
//   sprintf(buf, " %d (%.1f)\n", rep+1, timetoend);
//   Rcout <<  buf;
// }
//
// //[[Rcpp::export]]
// void endMcmcTimer() {
//   time_t ctime = time(NULL);
//   char buf[32];
//
//   sprintf(buf, " Total Time Elapsed: %.2f \n", difftime(ctime, itime) / 60.0);
//   Rcout << buf;
//
//   itime = 0;
// }
