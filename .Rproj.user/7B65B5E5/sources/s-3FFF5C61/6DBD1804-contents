#include "BFMediate.h"




// utility functions
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
double lldstar(arma::vec const& dstar, arma::vec const& y, vec const& mu, double ssq_y_tilde){       //y is y*, mu (= z) is y, ssq_y_tilde is the variance of indicator latent variables
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

