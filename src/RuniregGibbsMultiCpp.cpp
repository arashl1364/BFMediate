#include "BFMediate.h"


//This code just performs single runiregGibbs sampling R number of times
List RuniregGibbsMultiCpp(arma::mat const& y, arma::mat const& M, arma::vec const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq,
                          arma::vec sigmasq, int R, int keep, int nprint, bool betafix, bool sigmafix, arma::mat betavalue, arma::vec sigmavalue) {


  int N = X.size() ;
  int nvar = 3;
  //collection of moments
  arma::mat btildedraw(R,nvar);
  arma::cube IRdraw(nvar,nvar,R);
  arma::mat Xreg(N,3);
  Xreg.col(0) = ones(N);
  Xreg.col(2) = X;

  arma::mat draws;
  for (int rep=0; rep<R; rep++){

    Xreg.col(1) = trans(M.row(rep));


    List out = runiregGibbs_rcpp_me(trans(y.row(rep)) , Xreg, betabar, A, nu, ssq,
                         sigmasq(rep), 1, keep, nprint, betafix, sigmafix, betavalue, sigmavalue);

    btildedraw.row(rep) = as<arma::mat>(out["mubeta"]);
    IRdraw.slice(rep) = as<arma::cube>(out["IR"]);

  }

  return List::create(
    //moments collection for computing Bayes factor a la Savage Dickey (multiple indicator)
    Named("mubeta") = btildedraw,
    Named("IR")  = IRdraw);

}

