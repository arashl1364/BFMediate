// [[Rcpp::depends(RcppArmadillo)]]
//if you include bayesm.h here, you do not need to define functions before calling them
#include <RcppArmadillo.h>
// #include "rordprobitGibbs_me.cpp"
// #include "rordprobitGibbs_me_M.cpp"
// #include "rordprobitGibbs_me_MY.cpp"
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//MAIN FUNCTION-1--------------------------------------------------------------------------------------
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


//MAIN FUNCTION-2--------------------------------------------------------------------------------------
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




time_t itime;
char buf[100];


void startMcmcTimer() {
  itime = time(NULL);
  Rcout << " MCMC Iteration (est time to end - min) \n";
}


void infoMcmcTimer(int rep, int R) {
  time_t ctime = time(NULL);
  char buf[32];

  double timetoend = difftime(ctime, itime) / 60.0 * (R - rep - 1) / (rep+1);
  sprintf(buf, " %d (%.1f)\n", rep+1, timetoend);
  Rcout <<  buf;
}


void endMcmcTimer() {
  time_t ctime = time(NULL);
  char buf[32];

  sprintf(buf, " Total Time Elapsed: %.2f \n", difftime(ctime, itime) / 60.0);
  Rcout << buf;

  itime = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////            MAIN FUNCTION                  /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List Mediation_Ordered_Multi_Merr_cpp(arma::mat const& X, arma::mat const& m_star, arma::mat const& y_star, int k_M, int k_Y, int M_ind, int Y_ind,     //data
                                      arma::mat const& A_M, arma::vec const& betabar, arma::mat const& Ad_M, double s_M, arma::mat const& inc_root_M, arma::vec const& dstarbar_M, arma::vec const& betahat,               //priors_M
                                      arma::mat const& A_Y, arma::vec const& beta_2_bar, arma::mat const& Ad_Y, double s_Y, arma::mat const& inc_root_Y, arma::vec const& dstarbar_Y, arma::vec const& beta_2_hat,         //priors_Y
                                      int R, int keep, int nprint){
                                      // mat const& cutoff_M_init, mat const& M_tilde_init, vec const& beta_m_tilde_init, vec const& ssq_m_tilde_init, vec const& beta_init, vec const& M_init,                 //inital values_M
                                      // mat const& cutoff_Y_init, mat const& Y_tilde_init, vec const& beta_y_tilde_init, vec const& ssq_y_tilde_init, vec const& beta_2_init, vec const& Y_init){            //inital values_Y


  int mkeep;
  int ny = y_star.n_rows;

  int nvar_M = X.n_cols;
  int ncuts_M = k_M+1;
  // int ncut = ncuts-3;
  int ndstar_M = k_M-2;

  arma::mat betadraw(R/keep, nvar_M);
  arma::cube cutdraw_M(M_ind, ncuts_M, R/keep);
  arma::cube dstardraw_M(M_ind, ndstar_M,R/keep);
  arma::mat ssq_m_tilde_draw(R/keep, M_ind);
  arma::cube beta_m_tilde_draw(M_ind, 2, R/keep);
  arma::mat Mdraw(R/keep, ny);


  int nvar_Y = X.n_cols+1;
  int ncuts_Y = k_Y+1;
  // int ncut = ncuts-3;
  int ndstar_Y = k_Y-2;

  arma::mat beta_2_draw(R/keep, nvar_Y);
  arma::cube cutdraw_Y(Y_ind, ncuts_Y, R/keep);
  arma::cube dstardraw_Y(Y_ind, ndstar_Y,R/keep);
  arma::mat ssq_y_tilde_draw(R/keep, Y_ind);
  arma::cube beta_y_tilde_draw(Y_ind, 2, R/keep);
  arma::mat Ydraw(R/keep, ny);


  arma::mat mubeta_2_draw(R/keep, nvar_Y);
  arma::cube varbeta_2_draw(nvar_Y, nvar_Y, R/keep);


  // set Initial values
  arma::vec M = randu<vec>(ny);
  arma::mat XM(ny,nvar_Y);
  XM.col(0).ones();
  XM.col(1) = M;
  // XM.col(1) = M_init;       ///////// !!! CHANGE here after testing !!! ///////////
  XM.col(2) = X.col(1);

  arma::mat olddstar_M(M_ind,ndstar_M);
  olddstar_M.zeros();
  arma::vec oldbeta = betahat;
  // vec oldbeta = beta_init;        ///////// !!! CHANGE here after testing !!! ///////////
  arma::vec oldz_M = M;
  // vec oldz_M = M_init;         ///////// !!! CHANGE here after testing !!! ///////////
  arma::mat old_m_tilde(ny, M_ind);
  old_m_tilde.randn();
  // old_m_tilde = M_tilde_init;        ///////// !!! CHANGE here after testing !!! ///////////
  arma::vec old_ssq_m_tilde(M_ind);
  old_ssq_m_tilde.ones();
  // old_ssq_m_tilde = ssq_m_tilde_init;     // CHANGE HERE AFTER TEST
  arma::mat old_beta_m_tilde(M_ind,2);
  old_beta_m_tilde.col(1).ones();
  // old_beta_m_tilde.col(0) = beta_m_tilde_init;      // CHANGE HERE AFTER TEST


  arma::mat olddstar_Y(Y_ind,ndstar_Y);
  olddstar_Y.zeros();
  // olddstar_Y = dstar_Y_init;
  arma::vec oldbeta_2 = beta_2_hat;
  // vec oldbeta_2 = beta_2_init;   ///////// !!! CHANGE here after testing !!! ///////////
  arma::vec oldz_Y = randu<vec>(ny);
  // vec oldz_Y = Y_init;              ///////// !!! CHANGE here after testing !!! ///////////
  arma::mat old_y_tilde(ny, Y_ind);
  old_y_tilde.randn();
  // old_y_tilde = Y_tilde_init;        ///////// !!! CHANGE here after testing !!! ///////////
  vec old_ssq_y_tilde(Y_ind);
  old_ssq_y_tilde.ones();
  // old_ssq_y_tilde = ssq_y_tilde_init;     // CHANGE HERE AFTER TEST
  arma::mat old_beta_y_tilde(Y_ind,2);
  old_beta_y_tilde.col(1).ones();
  // old_beta_y_tilde.col(0) = beta_y_tilde_init;      // CHANGE HERE AFTER TEST

  for(int rep=0; rep<R; rep++){


    List out_Y = rordprobitGibbs_me(y_star, XM, k_Y, A_Y, beta_2_bar, Ad_Y,
                                    s_Y, inc_root_Y, dstarbar_Y, beta_2_hat, Y_ind,
                                    1, 1, 1,
                                    olddstar_Y, old_y_tilde, old_beta_y_tilde, old_ssq_y_tilde, oldbeta_2, oldz_Y); //, cutoff_Y_init);

    List out_M = rordprobitGibbs_me_M(oldz_Y, oldbeta_2, m_star, X, k_M, A_M, betabar, Ad_M,
                                      s_M, inc_root_M, dstarbar_M, betahat,M_ind,
                                      1, 1, 1,
                                      olddstar_M, old_m_tilde, old_beta_m_tilde, old_ssq_m_tilde, oldbeta, oldz_M); //, cutoff_M_init);

    //updating parameters

    //////// FIXING parameters //////////
    XM.col(1) = as<arma::vec>(out_M["zdraw"]);
    oldz_M = as<arma::vec>(out_M["zdraw"]);
    oldbeta = as<arma::vec>(out_M["betadraw"]);
    olddstar_M = as<arma::mat>(out_M["dstardraw"]);
    old_m_tilde = as<arma::mat>(out_M["y_tilde_draw"]);
    old_beta_m_tilde = as<arma::mat>(out_M["beta_tilde_draw"]);
    old_ssq_m_tilde = as<arma::vec>(out_M["ssq_y_tilde_draw"]);
    //
    oldz_Y = as<arma::vec>(out_Y["zdraw"]);
    oldbeta_2 = as<arma::vec>(out_Y["betadraw"]);
    olddstar_Y = as<arma::mat>(out_Y["dstardraw"]);
    old_y_tilde = as<arma::mat>(out_Y["y_tilde_draw"]);
    old_beta_y_tilde = as<arma::mat>(out_Y["beta_tilde_draw"]);
    old_ssq_y_tilde = as<arma::vec>(out_Y["ssq_y_tilde_draw"]);
    /////////////////////////////////////////////

    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      cutdraw_M.slice(mkeep-1) = as<arma::mat>(out_M["cutdraw"]);
      dstardraw_M.slice(mkeep-1) = as<arma::mat>(out_M["dstardraw"]);
      betadraw(mkeep-1,span::all) = trans(as<arma::vec>(out_M["betadraw"]));
      beta_m_tilde_draw.slice(mkeep-1) = as<arma::mat>(out_M["beta_tilde_draw"]);
      ssq_m_tilde_draw(mkeep-1,span::all) = trans(as<arma::vec>(out_M["ssq_y_tilde_draw"]));
      Mdraw(mkeep-1,span::all) = trans(as<arma::vec>(out_M["zdraw"]));

      cutdraw_Y.slice(mkeep-1) = as<arma::mat>(out_Y["cutdraw"]);
      dstardraw_Y.slice(mkeep-1) = as<arma::mat>(out_Y["dstardraw"]);
      beta_2_draw(mkeep-1,span::all) = trans(as<arma::vec>(out_Y["betadraw"]));
      beta_y_tilde_draw.slice(mkeep-1) = as<arma::mat>(out_Y["beta_tilde_draw"]);
      ssq_y_tilde_draw(mkeep-1,span::all) = trans(as<arma::vec>(out_Y["ssq_y_tilde_draw"]));
      Ydraw(mkeep-1,span::all) = trans(as<arma::vec>(out_Y["zdraw"]));

      mubeta_2_draw(mkeep-1,span::all) = trans(as<arma::vec>(out_Y["mubeta"]));
      varbeta_2_draw.slice(mkeep-1) = as<arma::mat>(out_Y["varbeta"]);
    }

  }

  return List::create(
    Named("cutdraw_M") = cutdraw_M,
    Named("dstardraw_M") = dstardraw_M,
    Named("betadraw") = betadraw,
    Named("beta_m_tilde_draw") = beta_m_tilde_draw,
    Named("ssq_m_tilde_draw") = ssq_m_tilde_draw,
    Named("Mdraw") = Mdraw,
    Named("cutdraw_Y") = cutdraw_Y,
    Named("dstardraw_Y") = dstardraw_Y,
    Named("beta_2_draw") = beta_2_draw,
    Named("beta_y_tilde_draw") = beta_y_tilde_draw,
    Named("ssq_y_tilde_draw") = ssq_y_tilde_draw,
    Named("Ydraw") = Ydraw,
    Named("mubeta_2_draw") = mubeta_2_draw,
    Named("varbeta_2_draw") = varbeta_2_draw
  );


}
