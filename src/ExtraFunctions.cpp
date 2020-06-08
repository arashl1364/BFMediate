#include "BFMediate.h"



// Function used in MeasurementMYCat.cpp
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


//in YSampler X is the matrix of (constant, X, M), y is the observed y_tilde, z is the latent Y and y_tilde is the latent y_star
List YSampler(arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
                        double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                        int const& Y_ind,
                        int R, int keep, int nprint,
                        arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz){

  int i;

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

  //beta moments for bayes factor
  arma::vec mubeta(nvar);
  arma::mat varbeta(nvar,nvar);

  //start main iteration loop
  for (int rep=0; rep<R; rep++){

    //draw gammas and y tilde's
    for(int ind=0; ind<Y_ind; ind++){

      //draw gamma
      metropout = dstarRwMetrop3(y(span::all,ind),beta_tilde(ind,0)+z,trans(olddstar(ind,span::all)),s,inc_root,dstarbar, rootdi, ncut, ssq_y_tilde[ind]);
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

////////////////////////////////////////////////////////////////////////////////////
//in YSampler X is the matrix of (constant, X), y is the observed m_tilde, z is the latent M and y_tilde is the latent m_star
List MSampler(arma::vec const& dep, arma::vec const& beta_2,  arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad,
                          double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat,
                          int const& Y_ind,
                          int R, int keep, int nprint,
                          arma::mat olddstar, arma::mat const& old_y_tilde, arma::mat const& old_beta_tilde, arma::vec const& old_ssq_y_tilde, arma::vec const& oldbeta, arma::vec const& oldz){

  int i;

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

  // start main iteration loop
  for (int rep=0; rep<R; rep++){


    //draw beta given z and rest             THIS REMAINS THE SAME P(beta_1|M,X)
    beta = breg1(root,X,z,Abetabar);

    for(int ind=0; ind<Y_ind; ind++){

      //draw gamma
      metropout = dstarRwMetrop3(y(span::all,ind),beta_tilde(ind,0)+z,trans(olddstar(ind,span::all)),s,inc_root,dstarbar, rootdi, ncut, ssq_y_tilde[ind]);
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
