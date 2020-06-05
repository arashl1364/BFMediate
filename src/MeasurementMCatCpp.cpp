#include "BFMediate.h"
// // [[Rcpp::depends(RcppArmadillo)]]
// // #include <bayesm.h> //if you include bayesm.h here, you do not need to define functions before calling them
// #include <RcppArmadillo.h>
// using namespace arma; // use the Armadillo library for matrix computations
// using namespace Rcpp;

// List runiregG(vec const& y, mat const& X, mat const& XpX, vec const& Xpy, double sigmasq, mat const& A,
//                 vec const& Abetabar, double nu, double ssq) {
//
//   // Keunwoo Kim 09/16/2014
//
//   // Purpose:
//   //  perform one Gibbs iteration for Univ Regression Model
//   //  only does one iteration so can be used in rhierLinearModel
//
//   // Model:
//   //  y = Xbeta + e  e ~N(0,sigmasq)
//   //  y is n x 1
//   //  X is n x k
//   //  beta is k x 1 vector of coefficients
//
//   // Prior:
//   //  beta ~ N(betabar,A^-1)
//   //  sigmasq ~ (nu*ssq)/chisq_nu
//
//   unireg out_struct;
//
//   int n = y.size();
//   int k = XpX.n_cols;
//
//   //first draw beta | sigmasq
//   mat IR = solve(trimatu(chol(XpX/sigmasq+A)), eye(k,k)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
//   vec btilde = (IR*trans(IR)) * (Xpy/sigmasq + Abetabar);
//   vec beta = btilde + IR*vec(rnorm(k));
//
//   //now draw sigmasq | beta
//   double s = sum(square(y-X*beta));
//   sigmasq = (s + nu*ssq)/rchisq(1,nu+n)[0]; //rchisq returns a vectorized object, so using [0] allows for the conversion to double
//
//   return List::create(
//     Named("beta") = beta,
//     Named("sigmasq") = sigmasq,
//     Named("mubeta") = btilde,
//     Named("varbeta") = IR*trans(IR));
// }



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//MAIN FUNCTION---------------------------------------------------------------------------------------
List MeasurementMCatCpp(arma::vec const& dep,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad, arma::mat const& A_2, arma::vec const& betabar_2,
                        double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                        int const& M_ind,
                        int R, int keep, int nprint){
  // mat const& cutoff_Y_init, mat const& Y_tilde_init, vec const& beta_tilde_init, vec const& ssq_y_tilde_init, vec const& beta_init, vec const& beta_2_init, vec const& Y_init){
  // vec const& z_init){

  // Modified by Arash Laghaie 16/10/2018

  // "dep" is the dependent variable in the mediation Y
  //  z is the continuous latent mediator M
  //  y_tilde is the continuous latent mediator with measurement error m_tilde  (REPLACED)
  //  y is the discrete mediator m*
  //  beta is beta_1 and beta_2 is beta_2


  // Keunwoo Kim 09/09/2014

  // Purpose: draw from posterior for ordered probit using Gibbs Sampler and metropolis RW

  // Arguments:
  //  Data
  //    X is nobs x nvar, y is nobs vector of 1,2,.,k (ordinal variable)
  //  Prior
  //    A is nvar x nvar prior preci matrix
  //    betabar is nvar x 1 prior mean
  //    Ad is ndstar x ndstar prior preci matrix of dstar (ncut is number of cut-offs being estimated)
  //    dstarbar is ndstar x 1 prior mean of dstar
  //  Mcmc
  //    R is number of draws
  //    keep is thinning parameter
  //    nprint - prints the estimated time remaining for every nprint'th draw
  //    s is scale parameter of random work Metropolis

  // Output: list of betadraws and cutdraws

  // Model:
  //    z=Xbeta + e  < 0  e ~N(0,1)
  //    y=1,..,k, if z~c(c[k], c[k+1])

  //    cutoffs = c[1],..,c[k+1]
  //    dstar = dstar[1],dstar[k-2]
  //    set c[1]=-100, c[2]=0, ...,c[k+1]=100

  //    c[3]=exp(dstar[1]),c[4]=c[3]+exp(dstar[2]),...,
  //    c[k]=c[k-1]+exp(datsr[k-2])

  // Note: 1. length of dstar = length of cutoffs - 3
  //       2. Be careful in assessing prior parameter, Ad.  .1 is too small for many applications.

  // Prior:
  //  beta ~ N(betabar,A^-1)
  //  beta_2 ~ N(beta_2_bar, A_2^-1)
  //  dstar ~ N(dstarbar, Ad^-1)

  // int stay;
  int i, mkeep;

  List metropout;

  int nvar = X.n_cols;     //X has 2 columns (first column is iota)
  int ncuts = k+1;
  int ncut = ncuts-3;
  int ndstar = k-2;
  int ny = y.n_rows;
  arma::vec z = zeros(ny);
  // vec z = Y_init;       // CHANGE HERE AFTER TEST

  arma::mat zdraw(R/keep,ny);
  arma::mat betadraw(R/keep, nvar);
  arma::mat beta_2_draw(R/keep, nvar+1);
  arma::mat ssq_m_tilde_draw(R/keep, M_ind);
  arma::cube lambdadraw(M_ind, 2, R/keep);
  arma::cube cutdraw(M_ind, ncuts, R/keep);
  arma::cube dstardraw(M_ind, ndstar,R/keep);
  arma::vec cutoff1(ny);
  arma::vec cutoff2(ny);
  arma::vec sigma(X.n_rows); sigma.ones();

  arma::vec mubeta_2_draw(R); //(R/keep,nvar+1);
  arma::vec varbeta_2_draw(R); //(nvar+1,nvar+1,R/keep);
  arma::vec ssq_Y_draw(R/keep);

  // compute the inverse of trans(X)*X+A
  arma::mat ucholinv = solve(trimatu(chol(trans(X)*X+A)), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
  // printf("442");
  arma::mat XXAinv = ucholinv*trans(ucholinv);

  arma::mat root = chol(XXAinv);
  // printf("446");
  arma::vec Abetabar = trans(A)*betabar;

  // compute the inverse of Ad
  ucholinv = solve(trimatu(chol(Ad)), eye(ndstar,ndstar));
  // printf("451");
  arma::mat Adinv = ucholinv*trans(ucholinv);

  arma::mat rootdi = chol(Adinv);
  // printf("455");
  // set initial values for MCMC
  mat olddstar(M_ind,ndstar);
  olddstar.zeros();
  // olddstar = dstar_init;
  arma::vec beta = betahat;
  // vec beta = beta_init;             //CHANGE HERE AFTER TEST
  arma::vec beta_2(3);
  // vec beta_2 = beta_2_init;         //CHANGE HERE AFTER TEST
  arma::mat cutoffs(M_ind,ncuts);
  arma::mat m_tilde(ny,M_ind);
  m_tilde.randn();
  // y_tilde = Y_tilde_init;     // CHANGE HERE AFTER TEST
  arma::vec ssq_m_tilde(M_ind);
  ssq_m_tilde.ones();
  // ssq_m_tilde = ssq_y_tilde_init;     // CHANGE HERE AFTER TEST
  arma::mat iota(ny,1);
  iota.ones();
  arma::mat iota_z(ny,2);
  iota_z.col(0) = iota.col(0);
  iota_z.col(1) = z;
  arma::mat lambda(M_ind,2);
  lambda.col(1).ones();    //coefficients of M in indicator equations are fixed to 1
  // lambda.col(0) = beta_tilde_init;      // CHANGE HERE AFTER TEST
  arma::vec betabar_tilde(1);  // we only need this for the regression when estimating the intercept of m_tilde
  betabar_tilde.zeros();
  arma::mat A_tilde(1,1);
  A_tilde(0,0) = .01;
  double ssq_Y = 1;

  for(int ind=0; ind<M_ind; ind++){

    cutoffs(ind,span::all) = trans(dstartoc(trans(olddstar(ind,span::all))));

  }
  // cutoffs = cutoff_Y_init;      // CHANGE HERE AFTER TEST

  //moments
  arma::vec mubeta(nvar+1);
  arma::mat varbeta(nvar+1,nvar+1);

  arma::mat XM(ny,nvar+1);
  XM.col(0).ones();
  XM.col(1) = z;
  XM.col(2) = X.col(1);

  // int first = 0;

  // start main iteration loop
  for (int rep=0; rep<R; rep++){

    // ////////////////////////////////////////////////////////////////
    // //draw beta_2 given Y, M, and X    (Conjugate regression with ssq_y=1)
    // XM.col(1) = z;
    // // compute the inverse of trans(X)*X+A
    // arma::mat ucholinv_2 = solve(trimatu(chol(trans(XM)*XM+A_2)), eye(nvar+1,nvar+1)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
    // // printf("505");
    // arma::mat XXAinv_2 = ucholinv_2*trans(ucholinv_2);
    //
    // arma::mat root_2 = chol(XXAinv_2);
    // // printf("509");
    // arma::vec Abetabar_2 = trans(A_2)*betabar_2;
    //
    // List beta_out = breg2(root_2,XM,dep,Abetabar_2);
    // beta_2 = as<vec>(beta_out["beta"]);
    // mubeta = as<vec>(beta_out["mubeta"]);
    // varbeta = as<mat>(beta_out["varbeta"]);
    //
    // arma::vec dep_tilde = dep - beta_2[0] - beta_2[2]*X.col(1);
    // ////////////////////////////////////////////////////////////////
    //draw beta_2 given Y, M, and X    (Gibbs step)
    XM.col(1) = z;
    arma::mat XpX = trans(XM)*XM;
    arma::vec Xpy = trans(XM)*dep;
    arma::vec Abetabar_2 = trans(A_2)*betabar_2;

    List beta_out = runiregG(dep, XM, XpX, Xpy, ssq_Y, A_2, Abetabar_2, 3, 1);   //last 2 arguments are nu and ssq

    beta_2 = as<vec>(beta_out["beta"]);
    ssq_Y =  as<double>(beta_out["sigmasq"]);
    mubeta = as<vec>(beta_out["mubeta"]);
    varbeta = as<mat>(beta_out["varbeta"]);
    ////////////////////////////////////////////////////////////////
    //draw beta given z and rest             THIS REMAINS THE SAME P(beta_1|M,X)
    beta = breg1(root,X,z,Abetabar);

    // if (first == 0) cutoffs = cutoffs_init;

    for(int ind=0; ind<M_ind; ind++){

      //draw gamma
      metropout = dstarRwMetrop(y(span::all,ind),lambda(ind,0)+z,trans(olddstar(ind,span::all)),s,inc_root,dstarbar, rootdi, ncut, ssq_m_tilde[ind]);
      olddstar(ind,span::all) = trans(as<vec>(metropout["dstardraw"])); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      cutoffs(ind,span::all) = trans(dstartoc(trans(olddstar(ind,span::all))));

      //draw m_tilde's
      arma::vec cutoff1_tilde(ny);
      arma::vec cutoff2_tilde(ny);
      arma::vec temp_sigma_tilde(ny);
      temp_sigma_tilde.fill(sqrt(ssq_m_tilde[ind]));
      for (i=0; i<ny; i++){
        cutoff1_tilde[i] = cutoffs(ind,y(i,ind)-1);  //lower bounds
        cutoff2_tilde[i] = cutoffs(ind,y(i,ind));    //upper bounds
      }
      m_tilde(span::all,ind) = rtrunVec(lambda(ind,0) + z, temp_sigma_tilde, cutoff1_tilde, cutoff2_tilde);

      //draw ssq_m_tilde and lambda (intercepts)
      if(ind == 0){
        iota_z.col(1) = z;
        List tilde_out = runiregGibbs_betafix(m_tilde(span::all,ind), iota_z, betabar_tilde, A_tilde, 3, 1, ssq_m_tilde[ind], 1, 1, 1, 1);//y, X, betabr, A, , nu, ssq, sigmasq, R, keep, nprint, betafix (we fix beta here)
        lambda(ind, 0) = 0;  //intercept of the first indicator m_tilde is fixed to 0
        ssq_m_tilde[ind] =  as<double>(tilde_out["sigmasqdraw"]);
      }
      else{
        List tilde_out = runiregGibbs_betafix(m_tilde(span::all,ind)-z, iota, betabar_tilde, A_tilde, 3, 1, ssq_m_tilde[ind], 1, 1, 1, 0);//y, X, betabr, A, , nu, ssq, sigmasq, R, keep, nprint, betafix (we don't fix beta here)
        lambda(ind, 0) = as<double>(tilde_out["betadraw"]);
        ssq_m_tilde[ind] =  as<double>(tilde_out["sigmasqdraw"]);
      }
    }

    //draw z given beta, beta_2, dep, cutoffs, y
    arma::vec p(M_ind+1);
    arma::mat q(M_ind+1,1);
    q(0,0) = beta_2(1)/sqrt(ssq_Y);    //in MeasurementMYCat  ssq_Y=1
    q(span(1,M_ind),0) = 1/sqrt(ssq_m_tilde);
    // compute the inverse of trans(X)*X+A where X is q, A is 1, and betabar is (beta_0 + M*beta_2 + X*beta_3)=X*beta
    arma::mat ucholinv_tilde = solve(trimatu(chol(trans(q)*q+1)), eye(1,1)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
    arma::mat XXAinv_tilde = ucholinv_tilde*trans(ucholinv_tilde);
    arma::mat root_tilde = chol(XXAinv_tilde);
    arma::vec Abetabar_tilde(1);
    for (i=0; i<ny; i++){
      p(0) = (dep(i) - beta_2(0) - beta_2(2)*X(i,1))/sqrt(ssq_Y);    //Y_i - beta_0 - beta_3*X_i    (likelihood1: Y|M,beta_2,beta_3,ssq_Y) , in MeasurementMYCat  ssq_Y=1
      p(span(1,M_ind)) = (trans(m_tilde.row(i)) - lambda.col(0))/sqrt(ssq_m_tilde);  //(likelihood2: m_tilde|M,lambda,ssq_m_tilde)
      Abetabar_tilde = X.row(i)*beta;  //A=1
      z(i) = conv_to<double>::from(breg1(root_tilde,q,p,Abetabar_tilde));

    }


    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      zdraw(mkeep-1,span::all) = trans(z);
      cutdraw.slice(mkeep-1) = cutoffs;
      dstardraw.slice(mkeep-1) = olddstar;
      betadraw(mkeep-1,span::all) = trans(beta);
      beta_2_draw(mkeep-1,span::all) = trans(beta_2);
      lambdadraw.slice(mkeep-1) = lambda;
      ssq_m_tilde_draw(mkeep-1,span::all) = trans(ssq_m_tilde);
      mubeta_2_draw(mkeep-1) = mubeta(nvar);  //(mkeep-1,span::all) = trans(mubeta);
      varbeta_2_draw(mkeep-1) = varbeta(nvar,nvar); //.slice(mkeep-1) = varbeta ;
      ssq_Y_draw(mkeep-1) = ssq_Y;
    }
  }
  // double accept = 1-sum(staydraw)/(R/keep);
  // if (nprint>0) endMcmcTimer();

  return List::create(
    Named("M") = zdraw,
    Named("cutoff_M") = cutdraw,
    // Named("dstardraw") = dstardraw,
    Named("beta_1") = betadraw,
    Named("beta_2") = beta_2_draw,
    Named("lambda") = lambdadraw,
    Named("ssq_m_star") = ssq_m_tilde_draw,
    Named("mu_draw") = mubeta_2_draw,
    Named("var_draw") = varbeta_2_draw,
    Named("ssq_Y") = ssq_Y_draw
  );
  // return List::create(
  //   Named("M_draw") = zdraw,
  //   Named("cutdraw") = cutdraw,
  //   Named("dstardraw") = dstardraw,
  //   Named("beta_1_draw") = betadraw,
  //   Named("beta_2_draw") = beta_2_draw,
  //   Named("beta_tilde_draw") = beta_tilde_draw,
  //   Named("ssq_y_tilde_draw") = ssq_y_tilde_draw,
  //   Named("mubeta_2_draw") = mubeta_2_draw,
  //   Named("varbeta_2_draw") = varbeta_2_draw,
  //   Named("ssq_Y_draw") = ssq_Y_draw
  // );
}






