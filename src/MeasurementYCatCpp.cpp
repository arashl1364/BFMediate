#include "BFMediate.h"

//MAIN FUNCTION---------------------------------------------------------------------------------------
List MeasurementYCatCpp(arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
                                            double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                                            int const& Y_ind,
                                            int R, int keep, int nprint){
                                       // mat const& cutoff_Y_init, mat const& Y_tilde_init, vec const& beta_tilde_init, vec const& ssq_y_tilde_init, vec const& beta_2_init, vec const& Y_init){

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
  //  dstar ~ N(dstarbar, Ad^-1)

  // int stay;
  int i, mkeep;

  List metropout;

  int nvar = X.n_cols;   //X has 3 columns (1,M,X)
  int ncuts = k+1;
  int ncut = ncuts-3;
  int ndstar = k-2;
  int ny = y.n_rows;
  arma::vec z = zeros(ny);
  // vec z = Y_init;       // CHANGE HERE AFTER TEST

  arma::mat zdraw(R/keep,ny);
  arma::mat betadraw(R/keep, nvar);
  arma::mat ssq_y_tilde_draw(R/keep, Y_ind);
  arma::cube taudraw(Y_ind, 2, R/keep);
  arma::cube cutdraw(Y_ind, ncuts, R/keep);
  arma::cube dstardraw(Y_ind, ndstar,R/keep);
  arma::vec cutoff1(ny);
  arma::vec cutoff2(ny);
  arma::vec sigma(X.n_rows); sigma.ones();

  arma::vec mubeta_2_draw(R);    //(R/keep,nvar);
  arma::vec varbeta_2_draw(R);   //(nvar,nvar,R/keep);

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
  arma::mat olddstar(Y_ind,ndstar);
  olddstar.zeros();
  arma::vec beta = betahat;
  arma::mat cutoffs(Y_ind,ncuts);
  arma::vec oldll(Y_ind);
  arma::mat y_tilde(ny,Y_ind);
  y_tilde.randn();
  // y_tilde = Y_tilde_init;     // CHANGE HERE AFTER TEST
  arma::vec ssq_y_tilde(Y_ind);
  ssq_y_tilde.ones();
  // ssq_y_tilde = ssq_y_tilde_init;     // CHANGE HERE AFTER TEST
  arma::mat iota(ny,1);
  iota.ones();
  arma::mat iota_z(ny,2);
  iota_z.col(0) = iota.col(0);
  iota_z.col(1) = z;
  arma::mat tau(Y_ind,2);
  tau.col(1).ones();
  // beta_tilde.col(0) = beta_tilde_init;      // CHANGE HERE AFTER TEST
  arma::vec betabar_tilde(1);  // we only need this for the regression when estimating the intercept of y_tilde
  betabar_tilde.zeros();
  arma::mat A_tilde(1,1);
  A_tilde(0,0) = .01;

  for(int ind=0; ind<Y_ind; ind++){

    cutoffs(ind,span::all) = trans(dstartoc(trans(olddstar(ind,span::all))));
    oldll[ind] = lldstar(trans(olddstar(ind,span::all)), y(span::all,ind), tau(ind,0)+z, ssq_y_tilde[ind]);

  }
  // cutoffs = cutoff_Y_init;      // CHANGE HERE AFTER TEST

  //beta moments for bayes factor
  arma::vec mubeta(nvar);
  arma::mat varbeta(nvar,nvar);

  // if (nprint>0) startMcmcTimer();

  //start main iteration loop
  for (int rep=0; rep<R; rep++){

    //draw gammas and y tilde's
    for(int ind=0; ind<Y_ind; ind++){

      //draw gamma
      metropout = dstarRwMetrop2(y(span::all,ind),tau(ind,0)+z,trans(olddstar(ind,span::all)),s,inc_root,dstarbar,oldll[ind],rootdi, ncut, ssq_y_tilde[ind]);
      olddstar(ind,span::all) = trans(as<arma::vec>(metropout["dstardraw"])); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
      oldll[ind] =  as<double>(metropout["oldll"]);
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
      y_tilde(span::all,ind) = rtrunVec(tau(ind,0) + z, temp_sigma_tilde, cutoff1_tilde, cutoff2_tilde);

      //draw ssq_y_tilde and tau (intercepts)
      if(ind == 0){
        iota_z.col(1) = z;
        List tilde_out = runiregGibbs_betafix(y_tilde(span::all,ind), iota_z, betabar_tilde, A_tilde, 3, 1, ssq_y_tilde[ind], 1, 1, 1, 1);//y, X, betabr, A, , nu, ssq, sigmasq, R, keep, nprint, betafix (we fix beta here)
        tau(ind, 0) = 0;  //intercept of the first indicator y_tilde is fixed to 0
        ssq_y_tilde[ind] =  as<double>(tilde_out["sigmasqdraw"]);
      }
      else{
        List tilde_out = runiregGibbs_betafix(y_tilde(span::all,ind)-z, iota, betabar_tilde, A_tilde, 3, 1, ssq_y_tilde[ind], 1, 1, 1, 0);//y, X, betabr, A, , nu, ssq, sigmasq, R, keep, nprint, betafix (we don't fix beta here)
        tau(ind, 0) = as<double>(tilde_out["betadraw"]);
        ssq_y_tilde[ind] =  as<double>(tilde_out["sigmasqdraw"]);
      }
    }

    //draw beta given z and rest
    List beta_out = breg2(root,X,z,Abetabar);
    beta = as<arma::vec>(beta_out["beta"]);
    mubeta = as<arma::vec>(beta_out["mubeta"]);
    varbeta = as<arma::mat>(beta_out["varbeta"]);

    //draw z given beta, sigma, y, cut-offs
    arma::vec p(Y_ind);
    arma::mat q(Y_ind,1);
    q.col(0) = 1/sqrt(ssq_y_tilde);
    // compute the inverse of trans(X)*X+A where X is q, A is 1, and betabar is (beta_0 + M*beta_2 + X*beta_3)=X*beta
    arma::mat ucholinv_tilde = solve(trimatu(chol(trans(q)*q+1)), eye(1,1)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
    arma::mat XXAinv_tilde = ucholinv_tilde*trans(ucholinv_tilde);
    arma::mat root_tilde = chol(XXAinv_tilde);
    arma::vec Abetabar_tilde(1);
    for (i=0; i<ny; i++){
      p = (trans(y_tilde.row(i)) - tau.col(0))/sqrt(ssq_y_tilde);
      Abetabar_tilde = X.row(i)*beta;  //A=1
      z(i) = conv_to<double>::from(breg1(root_tilde,q,p,Abetabar_tilde));

    }

    //print time to completion and draw # every nprint'th draw
    // if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);

    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      zdraw(mkeep-1,span::all) = trans(z);
      cutdraw.slice(mkeep-1) = cutoffs;
      dstardraw.slice(mkeep-1) = olddstar;
      betadraw(mkeep-1,span::all) = trans(beta);
      taudraw.slice(mkeep-1) = tau;
      ssq_y_tilde_draw(mkeep-1,span::all) = trans(ssq_y_tilde);
      mubeta_2_draw(mkeep-1) = mubeta(nvar-1);
      varbeta_2_draw(mkeep-1) = varbeta(nvar-1,nvar-1);

    }
  }
  // double accept = 1-sum(staydraw)/(R/keep);
  // if (nprint>0) endMcmcTimer();

  return List::create(
    Named("Y_draw") = zdraw,
    Named("cutoff_Y") = cutdraw,
    // Named("dstardraw") = dstardraw,
    Named("beta_Y") = betadraw,
    Named("tau") = taudraw,
    Named("ssq_y_tilde_draw") = ssq_y_tilde_draw,
    Named("mu_draw") = mubeta_2_draw,
    Named("var_draw") = varbeta_2_draw
  );
}







