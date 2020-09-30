#include "BFMediate.h"


// [[Rcpp::export]]
List rbprobitGibbs_rcpp_me(arma::vec const& y, arma::mat const& X, arma::vec const& Abetabar, arma::mat const& root,
                        arma::vec beta, arma::vec const& sigma, arma::vec const& trunpt, arma::vec const& above, int R, int keep, int nprint){

// Keunwoo Kim 09/09/2014
// Edited by A.Laghaie 2020

// Purpose: draw from posterior for binary probit using Gibbs Sampler

// Arguments:
//  X is nobs x nvar, y is nobs vector of 0,1
//  A is nvar x nvar prior preci matrix
//  betabar is nvar x 1 prior mean
//  R is number of draws
//  keep is thinning parameter
//  nprint - prints the estimated time remaining for every nprint'th draw

// Output: list of betadraws

// Model: y = 1 if  w=Xbeta+e>0  e~N(0,1)

// Prior: beta ~ N(betabar,A^-1)

  int mkeep;
  arma::vec mu;
  arma::vec z;
  arma::mat IR;
  arma::vec btilde;

  int nvar = X.n_cols;

  arma::mat betadraw(R/keep, nvar);
  arma::cube IRdraw(nvar,nvar,R/keep);
  arma::mat btildedraw(R/keep,nvar);
  // arma::mat zdraw(R/keep, X.n_rows);
  // arma::vec nudraw(R/keep);
  // arma::vec Sdraw(R/keep);
  arma::vec temp(nvar);


  //start main iteration loop
  for (int rep=0; rep<R; rep++){

    // draw z given beta(i-1)
    mu = X*beta;
    z = trunNorm_vec(mu, sigma, trunpt, above);
    beta = breg1(root, X, z, Abetabar);

    //store beta moments
    temp = (trans(root)*root)*(trans(X)*z+Abetabar);
    btildedraw.row(rep) = trans(temp);   //mean
    IRdraw.slice(rep) = trans(root)*root;  //covariance matrix


    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw(mkeep-1, span::all) = trans(beta);
    }
  }


  // return(trans((trans(root)*root)*(trans(X)*z+Abetabar)));
  return List::create(Named("betadraw") = betadraw,
                      Named("mu_beta") = btildedraw,
                      Named("IR")  = IRdraw);
}
