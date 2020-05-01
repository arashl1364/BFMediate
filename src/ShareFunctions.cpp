#include "BFMediate.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;

// These are the utility functions used in BFMediate

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




List breg_me(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A) {

  // Purpose: compute the posterior moments for linear regression, sigmasq=1.0

  // Output: draw from posterior

  // Model: y = Xbeta + e  e ~ N(0,I)

  // Prior:  beta ~ N(betabar,A^-1)

  int k = betabar.size();
  arma::mat RA = chol(A);
  // printf("98");
  arma::mat W = arma::join_cols(X, RA); //same as rbind(X,RA)     RA is U in Rossi slides
  arma::vec z = arma::join_cols(y, RA*betabar);  //z is the v matrix in Rossi slides
  arma::mat IR = solve(trimatu(chol(trans(W)*W)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  // printf("102");

  return List::create(
    Named("mubeta") = (IR*trans(IR))*(trans(W)*z),
    Named("IR")  = IR);
}



NumericVector rtrun(NumericVector const& mu, NumericVector const& sigma,
                    NumericVector const& a, NumericVector const& b){

  // Wayne Taylor 9/7/2014

  // function to draw from univariate truncated norm
  // a is vector of lower bounds for truncation
  // b is vector of upper bounds for truncation

  NumericVector FA = pnorm((a-mu)/sigma);
  NumericVector FB = pnorm((b-mu)/sigma);

  return(mu+sigma*qnorm(runif(mu.size())*(FB-FA)+FA));
}


//dstartoc is a fuction to transfer dstar to its cut-off value
arma::vec dstartoc(arma::vec const& dstar){
  int ndstar = dstar.size();
  arma::vec c(ndstar+3);
  c[0] = -100;
  c[1] = 0;
  c(span(2,ndstar+1)) = cumsum(exp(dstar));
  c[ndstar+2] = 100;

  return (c);
}


// compute conditional likelihood of data given cut-offs
double lldstar(arma::vec const& dstar, arma::vec const& y, arma::vec const& mu, double ssq_y_tilde){       //y is y*, mu (= z) is y, ssq_y_tilde is the variance of indicator latent variables
  arma::vec gamma = dstartoc(dstar);

  double sigma_y_tilde = sqrt(ssq_y_tilde);
  int ny = y.size();
  NumericVector gamma1(ny);
  NumericVector gamma2(ny);
  for (int i=0; i<ny; i++){
    gamma1[i] = gamma(y[i]);
    gamma2[i] = gamma(y[i]-1);
  }
  NumericVector temp = pnorm((gamma1-as<NumericVector>(wrap(mu)))/sigma_y_tilde)-pnorm((gamma2-as<NumericVector>(wrap(mu)))/sigma_y_tilde); //pnorm takes Rcpp type NumericVector, NOT arma objects of type vec
  arma::vec arg = as<arma::vec>(temp);
  double epsilon = 1.0/(10^-50);
  for (int j=0; j<ny; j++){
    if (arg[j]<epsilon){
      arg[j] = epsilon;
    }
  }
  return (sum(log(arg)));
}


double lndMvn(arma::vec const& x, arma::vec const& mu, arma::mat const& rooti){

  //Wayne Taylor 9/7/2014

  // function to evaluate log of MV Normal density with  mean mu, var Sigma
  // Sigma=t(root)%*%root   (root is upper tri cholesky root)
  // Sigma^-1=rooti%*%t(rooti)
  // rooti is in the inverse of upper triangular chol root of sigma
  //          note: this is the UL decomp of sigmai not LU!
  //                Sigma=root'root   root=inv(rooti)

  arma::vec z = vectorise(trans(rooti)*(x-mu));

  return((-(x.size()/2.0)*log(2*M_PI) -.5*(trans(z)*z) + sum(log(diagvec(rooti))))[0]);
}


List dstarRwMetrop(arma::vec const& y, arma::vec const& mu, arma::vec const& olddstar, double s, arma::mat const& inc_root,
                   arma::vec const& dstarbar, arma::mat const& rootdi, int ncut, double ssq_y_tilde){

  // function to execute rw metropolis for the dstar
  // y is n vector with element = 1,...,j
  // X is n x k matrix of x values
  // RW increments are N(0,s^2*t(inc.root)%*%inc.root)
  // prior on dstar is N(dstarbar,Sigma)  Sigma^-1=rootdi*t(rootdi)
  //  inc.root, rootdi are upper triangular
  //  this means that we are using the UL decomp of Sigma^-1 for prior
  // olddstar is the current
  //
  // int stay = 0;
  double unif;
  arma::vec dstardraw;

  arma::vec dstarc = olddstar + s*trans(inc_root)*vec(rnorm(ncut));
  double oldll = lldstar(olddstar, y, mu, ssq_y_tilde);
  double cll = lldstar(dstarc, y, mu, ssq_y_tilde);
  double clpost = cll + lndMvn(dstarc, dstarbar, rootdi);
  double ldiff = clpost - oldll - lndMvn(olddstar, dstarbar, rootdi);
  double alpha = exp(ldiff);

  if (alpha>1){
    alpha = 1.0;
  }

  if (alpha<1){
    unif = runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double by extracting the first element
  }
  else{
    unif = 0;
  }

  if (unif<=alpha){
    dstardraw = dstarc;
    oldll = cll;
  }
  else{
    dstardraw = olddstar;
    // stay = 1;
  }

  return List::create(
    Named("dstardraw") = dstardraw,
    Named("oldll") = oldll
    // Named("stay") = stay
  );
}



List dstarRwMetrop_1(arma::vec const& y, arma::vec const& mu, arma::vec const& olddstar, double s, arma::mat const& inc_root,
                     arma::vec const& dstarbar, arma::mat const& rootdi, int ncut, double ssq_y_tilde){

  // function to execute rw metropolis for the dstar
  // y is n vector with element = 1,...,j
  // X is n x k matrix of x values
  // RW increments are N(0,s^2*t(inc.root)%*%inc.root)
  // prior on dstar is N(dstarbar,Sigma)  Sigma^-1=rootdi*t(rootdi)
  //  inc.root, rootdi are upper triangular
  //  this means that we are using the UL decomp of Sigma^-1 for prior
  // olddstar is the current

  double unif;
  arma::vec dstardraw;

  arma::vec dstarc = olddstar + s*trans(inc_root)*vec(rnorm(ncut));
  double oldll = lldstar(olddstar, y, mu, ssq_y_tilde);
  double cll = lldstar(dstarc, y, mu, ssq_y_tilde);
  double clpost = cll + lndMvn(dstarc, dstarbar, rootdi);
  double ldiff = clpost - oldll - lndMvn(olddstar, dstarbar, rootdi);
  double alpha = exp(ldiff);

  if (alpha>1){
    alpha = 1.0;
  }

  if (alpha<1){
    unif = runif(1)[0]; //runif returns a NumericVector, so using [0] allows for conversion to double by extracting the first element
  }
  else{
    unif = 0;
  }

  if (unif<=alpha){
    dstardraw = dstarc;
  }
  else{
    dstardraw = olddstar;
  }

  return List::create(
    Named("dstardraw") = dstardraw
  );
}


List breg2(arma::mat const& root, arma::mat const& X, arma::vec const& y, arma::vec const& Abetabar) {

  // Arash Laghaie 10/20/2018
  // Equivalent to breg1 and additionally stores the beta moments

  arma::mat cov = trans(root)*root;

  return List::create(
    Named("beta") = cov*(trans(X)*y+Abetabar) + trans(root)*vec(rnorm(root.n_cols)),
    Named("mubeta") = cov*(trans(X)*y+Abetabar),
    Named("varbeta") = cov   //beta covariance matrix
  );
}





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
double norm_rs(double a, double b)
{
  double  x;
  x = Rf_rnorm(0.0, 1.0);
  while( (x < a) || (x > b) ) x = norm_rand();
  return x;
}

double half_norm_rs(double a, double b)
{
  double   x;
  x = fabs(norm_rand());
  while( (x<a) || (x>b) ) x = fabs(norm_rand());
  return x;
}

double unif_rs(double a, double b)
{
  double xstar, logphixstar, x, logu;

  // Find the argmax (b is always >= 0)
  // This works because we want to sample from N(0,1)
  if(a <= 0.0) xstar = 0.0;
  else xstar = a;
  logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);

  x = R::runif(a, b);
  logu = log(R::runif(0.0, 1.0));
  while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
  {
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
  }
  return x;
}

double exp_rs(double a, double b)
{
  double  z, u, rate;

  //  Rprintf("in exp_rs");
  rate = 1/a;
  //1/a

  // Generate a proposal on (0, b-a)
  z = R::rexp(rate);
  while(z > (b-a)) z = R::rexp(rate);
  u = R::runif(0.0, 1.0);

  while( log(u) > (-0.5*z*z))
  {
    z = R::rexp(rate);
    while(z > (b-a)) z = R::rexp(rate);
    u = R::runif(0.0,1.0);
  }
  return(z+a);
}


// NumericVector rtnm(Rcpp::NumericVector mus, Rcpp::NumericVector sigmas, Rcpp::NumericVector lower, Rcpp::NumericVector upper){
arma::vec rtnm(arma::vec mus, arma::vec sigmas, arma::vec lower, arma::vec upper){
  // omp_set_num_threads(cores);
  int nobs = mus.size();
  arma::vec out(nobs);
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double a,b, z, tmp, lograt;

  int  change;

  // #pragma omp parallel for schedule(dynamic)
  for(int i=0;i<nobs;i++) {

    a = (lower(i) - mus(i))/sigmas(i);
    b = (upper(i) - mus(i))/sigmas(i);
    change=0;
    // First scenario
    if( (a == R_NegInf) || (b == R_PosInf))
    {
      if(a == R_NegInf)
      {
        change = 1;
        a = -b;
        b = R_PosInf;
      }

      // The two possibilities for this scenario
      if(a <= 0.45) z = norm_rs(a, b);
      else z = exp_rs(a, b);
      if(change) z = -z;
    }
    // Second scenario
    else if((a * b) <= 0.0)
    {
      // The two possibilities for this scenario
      if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
      {
        z = norm_rs(a, b);
      }
      else z = unif_rs(a,b);
    }

    // Third scenario
    else
    {
      if(b < 0)
      {
        tmp = b; b = -a; a = -tmp; change = 1;
      }

      lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
      if(lograt <= logt2) z = unif_rs(a,b);
      else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
      else z = exp_rs(a,b);
      if(change) z = -z;
    }
    out(i)=sigmas(i)*z + mus(i);
  }

  return(out);
}






//The functions below are used to print the output from MCMC draws for many of the bayesm functions
arma::vec rtrunVec(arma::vec const& mu,arma::vec const& sigma, arma::vec const& a, arma::vec const& b){

  // Keunwoo Kim 06/20/2014

  //function to draw from univariate truncated norm
  //a is vector of lower bounds for truncation
  //b is vector of upper bounds for truncation

  int n = mu.size();
  arma::vec FA(n);
  arma::vec FB(n);
  arma::vec out(n);
  for (int i=0; i<n; i++) {
    FA[i] = R::pnorm((a[i]-mu[i])/sigma[i],0,1,1,0);
    FB[i] = R::pnorm((b[i]-mu[i])/sigma[i],0,1,1,0);
    out[i] = mu[i]+sigma[i]*R::qnorm(R::runif(0,1)*(FB[i]-FA[i])+FA[i],0,1,1,0);
  }

  return(out);
}



arma::vec breg1(arma::mat const& root, arma::mat const& X, arma::vec const& y, arma::vec const& Abetabar) {

  // Keunwoo Kim 06/20/2014

  // Purpose: draw from posterior for linear regression, sigmasq=1.0

  // Arguments:
  //  root = chol((X'X+A)^-1)
  //  Abetabar = A*betabar

  // Output: draw from posterior

  // Model: y = Xbeta + e  e ~ N(0,I)

  // Prior: beta ~ N(betabar,A^-1)

  arma::mat cov = trans(root)*root;

  return (cov*(trans(X)*y+Abetabar) + trans(root)*vec(rnorm(root.n_cols)));
}





//MAIN FUNCTION-1  For Mediation_Ordered_Multi_Merr.cpp--------------------------------------------------------------------------------------
List rordprobitGibbs_me(arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
                        double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                        int const& Y_ind,
                        int R, int keep, int nprint,
                        arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz){
  // mat const& cutoff_Y_init){




  int i; //j, mkeep;

  List metropout;

  int nvar = X.n_cols;
  int ncuts = k+1;
  int ncut = ncuts-3;
  int ndstar = k-2;
  int ny = y.n_rows;


  arma::mat zdraw(R/keep,ny);
  arma::mat betadraw(R/keep, nvar);
  arma::mat ssq_y_tilde_draw(R/keep, Y_ind);
  arma::cube beta_m_tilde_draw(Y_ind, 2, R/keep);
  arma::cube cutdraw(Y_ind, ncuts, R/keep);
  arma::cube dstardraw(Y_ind, ndstar,R/keep);
  arma::vec cutoff1(ny);
  arma::vec cutoff2(ny);
  arma::vec sigma(X.n_rows); sigma.ones();

  arma::mat mubetadraw(R/keep,nvar);
  arma::cube varbetadraw(nvar,nvar,R/keep);

  // compute the inverse of trans(X)*X+A
  arma::mat ucholinv = solve(trimatu(chol(trans(X)*X+A)), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  arma::mat XXAinv = ucholinv*trans(ucholinv);

  arma::mat root = chol(XXAinv);
  arma::vec Abetabar = trans(A)*betabar;

  // compute the inverse of Ad
  ucholinv = solve(trimatu(chol(Ad)), eye(ndstar,ndstar));
  arma::mat Adinv = ucholinv*trans(ucholinv);

  arma::mat rootdi = chol(Adinv);

  // set initial values for MCMC
  arma::vec beta = oldbeta;
  arma::mat cutoffs(Y_ind,ncuts);
  arma::mat y_tilde = old_y_tilde;
  arma::mat beta_tilde = old_beta_tilde;
  arma::vec ssq_y_tilde = old_ssq_y_tilde;
  arma::vec z = oldz;
  arma::mat iota(ny,1);
  iota.ones();
  arma::mat iota_z(ny,2);
  iota_z.col(0) = iota.col(0);
  iota_z.col(1) = z;

  arma::vec betabar_tilde(1);  // we only need this for the regression when estimating the intercept of y_tilde
  betabar_tilde.zeros();
  arma::mat A_tilde(1,1);
  A_tilde(0,0) = .01;

  for(int ind=0; ind<Y_ind; ind++) cutoffs(ind,span::all) = trans(dstartoc(trans(olddstar(ind,span::all))));

  // cutoffs = cutoff_Y_init;      // CHANGE HERE AFTER TEST

  //beta moments for bayes factor
  arma::vec mubeta(nvar);
  arma::mat varbeta(nvar,nvar);

  //start main iteration loop
  for (int rep=0; rep<R; rep++){

    //draw gammas and y tilde's
    for(int ind=0; ind<Y_ind; ind++){

      //draw gamma
      metropout = dstarRwMetrop_1(y(span::all,ind),beta_tilde(ind,0)+z,trans(olddstar(ind,span::all)),s,inc_root,dstarbar, rootdi, ncut, ssq_y_tilde[ind]);
      olddstar(ind,span::all) = trans(as<vec>(metropout["dstardraw"])); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      cutoffs(ind,span::all) = trans(dstartoc(trans(olddstar(ind,span::all))));

      //draw y_tilde's
      arma::vec cutoff1_tilde(ny);
      arma::vec cutoff2_tilde(ny);
      arma::vec temp_sigma_tilde(ny);
      temp_sigma_tilde.fill(sqrt(ssq_y_tilde[ind]));
      for (i=0; i<ny; i++){
        cutoff1_tilde[i] = cutoffs(ind,y(i,ind)-1);  //lower bounds
        cutoff2_tilde[i] = cutoffs(ind,y(i,ind));    //upper bounds
      }
      y_tilde(span::all,ind) = rtrunVec(beta_tilde(ind,0) + z, temp_sigma_tilde, cutoff1_tilde, cutoff2_tilde);

      //draw ssq_y_tilde and beta_tilde (intercepts)
      if(ind == 0){
        iota_z.col(1) = z;
        List tilde_out = runiregGibbs_betafix(y_tilde(span::all,ind), iota_z, betabar_tilde, A_tilde, 3, 1, ssq_y_tilde[ind], 1, 1, 1, 1);//y, X, betabr, A, , nu, ssq, sigmasq, R, keep, nprint, betafix (we fix beta here)
        beta_tilde(ind, 0) = 0;  //intercept of the first indicator y_tilde is fixed to 0
        ssq_y_tilde[ind] =  as<double>(tilde_out["sigmasqdraw"]);
      }
      else{
        List tilde_out = runiregGibbs_betafix(y_tilde(span::all,ind)-z, iota, betabar_tilde, A_tilde, 3, 1, ssq_y_tilde[ind], 1, 1, 1, 0);//y, X, betabr, A, , nu, ssq, sigmasq, R, keep, nprint, betafix (we don't fix beta here)
        beta_tilde(ind, 0) = as<double>(tilde_out["betadraw"]);
        ssq_y_tilde[ind] =  as<double>(tilde_out["sigmasqdraw"]);
      }
    }

    //draw beta given z and rest
    List beta_out = breg2(root,X,z,Abetabar);
    beta = as<arma::vec>(beta_out["beta"]);
    mubeta = as<arma::vec>(beta_out["mubeta"]);
    varbeta = as<arma::mat>(beta_out["varbeta"]);

    // draw z given beta, sigma, y, cut-offs
    arma::vec p(Y_ind);
    arma::mat q(Y_ind,1);
    q.col(0) = 1/sqrt(ssq_y_tilde);
    // compute the inverse of trans(X)*X+A where X is q, A is 1, and betabar is (beta_0 + M*beta_2 + X*beta_3)=X*beta
    arma::mat ucholinv_tilde = solve(trimatu(chol(trans(q)*q+1)), eye(1,1)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
    arma::mat XXAinv_tilde = ucholinv_tilde*trans(ucholinv_tilde);
    arma::mat root_tilde = chol(XXAinv_tilde);
    arma::vec Abetabar_tilde(1);
    for (i=0; i<ny; i++){
      p = (trans(y_tilde.row(i)) - beta_tilde.col(0))/sqrt(ssq_y_tilde);
      Abetabar_tilde = X.row(i)*beta;  //A=1
      z(i) = conv_to<double>::from(breg1(root_tilde,q,p,Abetabar_tilde));

    }
  }

  return List::create(
    Named("zdraw") = z,
    Named("cutdraw") = cutoffs,
    Named("dstardraw") = olddstar,
    Named("betadraw") = beta,
    Named("y_tilde_draw") = y_tilde,
    Named("beta_tilde_draw") = beta_tilde,
    Named("ssq_y_tilde_draw") = ssq_y_tilde,
    Named("mubeta") = mubeta,
    Named("varbeta") = varbeta
  );
}



//MAIN FUNCTION-2 For Mediation_Ordered_Multi_Merr.cpp-------------------------------------------------------------------------------------
List rordprobitGibbs_me_M(arma::vec const& dep, arma::vec const& beta_2,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
                          double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                          int const& Y_ind,
                          int R, int keep, int nprint,
                          arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz){
  // mat const& cutoff_M_init){


  int i;  //, mkeep;

  List metropout;

  int nvar = X.n_cols;
  int ncuts = k+1;
  int ncut = ncuts-3;
  int ndstar = k-2;
  int ny = y.n_rows;

  arma::mat zdraw(R/keep,ny);
  arma::mat betadraw(R/keep, nvar);
  arma::mat ssq_y_tilde_draw(R/keep, Y_ind);
  arma::cube beta_tilde_draw(Y_ind, 2, R/keep);
  arma::cube cutdraw(Y_ind, ncuts, R/keep);
  arma::cube dstardraw(Y_ind, ndstar,R/keep);
  arma::vec cutoff1(ny);
  arma::vec cutoff2(ny);
  arma::vec sigma(X.n_rows); sigma.ones();

  // compute the inverse of trans(X)*X+A
  arma::mat ucholinv = solve(trimatu(chol(trans(X)*X+A)), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  arma::mat XXAinv = ucholinv*trans(ucholinv);

  arma::mat root = chol(XXAinv);
  arma::vec Abetabar = trans(A)*betabar;

  // compute the inverse of Ad
  ucholinv = solve(trimatu(chol(Ad)), eye(ndstar,ndstar));
  arma::mat Adinv = ucholinv*trans(ucholinv);

  arma::mat rootdi = chol(Adinv);

  // set initial values for MCMC
  arma::vec beta = oldbeta;
  arma::mat cutoffs(Y_ind,ncuts);
  arma::mat y_tilde = old_y_tilde;
  arma::mat beta_tilde = old_beta_tilde;
  arma::vec ssq_y_tilde = old_ssq_y_tilde;
  arma::vec z = oldz;
  arma::mat iota(ny,1);
  iota.ones();
  arma::mat iota_z(ny,2);
  iota_z.col(0) = iota.col(0);
  iota_z.col(1) = z;

  arma::vec betabar_tilde(1);  // we only need this for the regression when estimating the intercept of y_tilde
  betabar_tilde.zeros();
  arma::mat A_tilde(1,1);
  A_tilde(0,0) = .01;

  for(int ind=0; ind<Y_ind; ind++) cutoffs(ind,span::all) = trans(dstartoc(trans(olddstar(ind,span::all))));

  // cutoffs = cutoff_M_init;      // CHANGE HERE AFTER TEST

  // start main iteration loop
  for (int rep=0; rep<R; rep++){


    //draw beta given z and rest             THIS REMAINS THE SAME P(beta_1|M,X)
    beta = breg1(root,X,z,Abetabar);

    for(int ind=0; ind<Y_ind; ind++){

      //draw gamma
      metropout = dstarRwMetrop_1(y(span::all,ind),beta_tilde(ind,0)+z,trans(olddstar(ind,span::all)),s,inc_root,dstarbar, rootdi, ncut, ssq_y_tilde[ind]);
      olddstar(ind,span::all) = trans(as<arma::vec>(metropout["dstardraw"])); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      cutoffs(ind,span::all) = trans(dstartoc(trans(olddstar(ind,span::all))));

      //draw y_tilde's
      arma::vec cutoff1_tilde(ny);
      arma::vec cutoff2_tilde(ny);
      arma::vec temp_sigma_tilde(ny);
      temp_sigma_tilde.fill(sqrt(ssq_y_tilde[ind]));
      for (i=0; i<ny; i++){
        cutoff1_tilde[i] = cutoffs(ind,y(i,ind)-1);  //lower bounds
        cutoff2_tilde[i] = cutoffs(ind,y(i,ind));    //upper bounds
      }
      y_tilde(span::all,ind) = rtrunVec(beta_tilde(ind,0) + z, temp_sigma_tilde, cutoff1_tilde, cutoff2_tilde);

      //draw ssq_y_tilde and beta_tilde (intercepts)
      if(ind == 0){
        iota_z.col(1) = z;
        List tilde_out = runiregGibbs_betafix(y_tilde(span::all,ind), iota_z, betabar_tilde, A_tilde, 3, 1, ssq_y_tilde[ind], 1, 1, 1, 1);//y, X, betabr, A, , nu, ssq, sigmasq, R, keep, nprint, betafix (we fix beta here)
        beta_tilde(ind, 0) = 0;  //intercept of the first indicator y_tilde is fixed to 0
        ssq_y_tilde[ind] =  as<double>(tilde_out["sigmasqdraw"]);
      }
      else{
        List tilde_out = runiregGibbs_betafix(y_tilde(span::all,ind)-z, iota, betabar_tilde, A_tilde, 3, 1, ssq_y_tilde[ind], 1, 1, 1, 0);//y, X, betabr, A, , nu, ssq, sigmasq, R, keep, nprint, betafix (we don't fix beta here)
        beta_tilde(ind, 0) = as<double>(tilde_out["betadraw"]);
        ssq_y_tilde[ind] =  as<double>(tilde_out["sigmasqdraw"]);
      }
    }

    //draw z given beta, beta_2, dep, cutoffs, y
    arma::vec p(Y_ind+1);
    arma::mat q(Y_ind+1,1);
    q(0,0) = beta_2(1);
    q(span(1,Y_ind),0) = 1/sqrt(ssq_y_tilde);
    // compute the inverse of trans(X)*X+A where X is q, A is 1, and betabar is (beta_0 + M*beta_2 + X*beta_3)=X*beta
    arma::mat ucholinv_tilde = solve(trimatu(chol(trans(q)*q+1)), eye(1,1)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
    arma::mat XXAinv_tilde = ucholinv_tilde*trans(ucholinv_tilde);
    arma::mat root_tilde = chol(XXAinv_tilde);
    arma::vec Abetabar_tilde(1);
    for (i=0; i<ny; i++){
      p(0) = dep(i) - beta_2(0) - beta_2(2)*X(i,1);    //Y_i - beta_0 - beta_3*X_i
      p(span(1,Y_ind)) = (trans(y_tilde.row(i)) - beta_tilde.col(0))/sqrt(ssq_y_tilde);
      Abetabar_tilde = X.row(i)*beta;  //A=1
      z(i) = conv_to<double>::from(breg1(root_tilde,q,p,Abetabar_tilde));

    }
  }

  return List::create(
    Named("zdraw") = z,
    Named("cutdraw") = cutoffs,
    Named("dstardraw") = olddstar,
    Named("betadraw") = beta,
    Named("y_tilde_draw") = y_tilde,
    Named("beta_tilde_draw") = beta_tilde,
    Named("ssq_y_tilde_draw") = ssq_y_tilde
  );
}

