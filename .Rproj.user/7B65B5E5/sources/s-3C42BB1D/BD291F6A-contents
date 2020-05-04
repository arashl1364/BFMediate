#include "BFMediate.h"



//EXTRA FUNCTIONS SPECIFIC TO THE MAIN FUNCTION--------------------------------------------
List runiregGibbs_betafix(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq,
                          double sigmasq, int R, int keep, int nprint, int betafix) {

  // Keunwoo Kim 09/09/2014

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
  arma::vec z, btilde,beta;

  int nvar = X.n_cols;
  int nobs = y.size();

  arma::vec sigmasqdraw(R/keep);
  arma::mat betadraw(R/keep, nvar);

  arma::mat XpX = trans(X)*X;
  arma::vec Xpy = trans(X)*y;

  arma::vec Abetabar = A*betabar;

  // if (nprint>0) startMcmcTimer();

  for (int rep=0; rep<R; rep++){

    if(betafix==1) beta <<0<<1;

    else{
      //first draw beta | sigmasq
      IR = solve(trimatu(chol(XpX/sigmasq+A)), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
      // printf("59");
      btilde = (IR*trans(IR)) * (Xpy/sigmasq+Abetabar);
      beta = btilde + IR*vec(rnorm(nvar));
    }

    //now draw sigmasq | beta
    s = sum(square(y-X*beta));
    sigmasq = (nu*ssq+s) / rchisq(1,nu+nobs)[0]; //rchisq returns a vectorized object, so using [0] allows for the conversion to double

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
    Named("sigmasqdraw") = NumericVector(sigmasqdraw.begin(),sigmasqdraw.end()));
}
