#include "BFMediate.h"


///////////////////////////////////////////////////////////////////
///////////////            MAIN FUNCTION          /////////////////
///////////////////////////////////////////////////////////////////
// ssq_m_star and ssq_y_star are called ssq_m_tilde and ssq_y_tilde
List MeasurementMYCatCpp(arma::mat const& X, arma::mat const& m_tilde, arma::mat const& y_tilde, int k_M, int k_Y, int M_ind, int Y_ind,     //data
                                      arma::mat const& A_M, arma::vec const& betabar, arma::mat const& Ad_M, double s_M, arma::mat const& inc_root_M, arma::vec const& dstarbar_M, arma::vec const& betahat,               //priors_M
                                      arma::mat const& A_Y, arma::vec const& beta_2_bar, arma::mat const& Ad_Y, double s_Y, arma::mat const& inc_root_Y, arma::vec const& dstarbar_Y, arma::vec const& beta_2_hat,         //priors_Y
                                      int R, int keep, int nprint){


  int mkeep;
  int ny = y_tilde.n_rows;

  int nvar_M = X.n_cols;
  int ncuts_M = k_M+1;
  int ndstar_M = k_M-2;

  arma::mat betadraw(R/keep, nvar_M);
  arma::cube cutoff_M(M_ind, ncuts_M, R/keep);
  arma::cube dstardraw_M(M_ind, ndstar_M,R/keep);
  arma::mat ssq_m_tilde_draw(R/keep, M_ind);
  arma::cube lambdadraw(M_ind, 2, R/keep);
  arma::mat Mdraw(R/keep, ny);


  int nvar_Y = X.n_cols+1;
  int ncuts_Y = k_Y+1;
  int ndstar_Y = k_Y-2;

  arma::mat beta_2_draw(R/keep, nvar_Y);
  arma::cube cutoff_Y(Y_ind, ncuts_Y, R/keep);
  arma::cube dstardraw_Y(Y_ind, ndstar_Y,R/keep);
  arma::mat ssq_y_tilde_draw(R/keep, Y_ind);
  arma::cube taudraw(Y_ind, 2, R/keep);
  arma::mat Ydraw(R/keep, ny);


  arma::vec mubeta_2_draw(R/keep);   //(R/keep, nvar_Y);
  arma::vec varbeta_2_draw(R/keep);  //(nvar_Y, nvar_Y, R/keep);


  // set Initial values
  arma::vec M = randu<arma::vec>(ny);
  arma::mat XM(ny,nvar_Y);
  XM.col(0).ones();
  XM.col(1) = M;
  XM.col(2) = X.col(1);

  arma::mat olddstar_M(M_ind,ndstar_M);
  olddstar_M.zeros();
  arma::vec oldbeta = betahat;
  arma::vec oldz_M = M;
  arma::mat old_m_tilde(ny, M_ind);
  old_m_tilde.randn();
  arma::vec old_ssq_m_tilde(M_ind);
  old_ssq_m_tilde.ones();
  arma::mat old_lambda(M_ind,2);
  old_lambda.col(1).ones();


  arma::mat olddstar_Y(Y_ind,ndstar_Y);
  olddstar_Y.zeros();
  arma::vec oldbeta_2 = beta_2_hat;
  arma::vec oldz_Y = randu<vec>(ny);
  arma::mat old_y_tilde(ny, Y_ind);
  old_y_tilde.randn();
  arma::vec old_ssq_y_tilde(Y_ind);
  old_ssq_y_tilde.ones();
  arma::mat old_tau(Y_ind,2);
  old_tau.col(1).ones();

  for(int rep=0; rep<R; rep++){


    List out_Y =  YSampler(y_tilde, XM, k_Y, A_Y, beta_2_bar, Ad_Y,
                                    s_Y, inc_root_Y, dstarbar_Y, beta_2_hat, Y_ind,
                                    1, 1, 1,
                                    olddstar_Y, old_y_tilde, old_tau, old_ssq_y_tilde, oldbeta_2, oldz_Y); //, cutoff_Y_init);

    List out_M =  MSampler(oldz_Y, oldbeta_2, m_tilde, X, k_M, A_M, betabar, Ad_M,
                                      s_M, inc_root_M, dstarbar_M, betahat,M_ind,
                                      1, 1, 1,
                                      olddstar_M, old_m_tilde, old_lambda, old_ssq_m_tilde, oldbeta, oldz_M); //, cutoff_M_init);

    //updating parameters
    XM.col(1) = as<arma::vec>(out_M["zdraw"]);
    oldz_M = as<arma::vec>(out_M["zdraw"]);
    oldbeta = as<arma::vec>(out_M["betadraw"]);
    olddstar_M = as<arma::mat>(out_M["dstardraw"]);
    old_m_tilde = as<arma::mat>(out_M["y_tilde_draw"]);
    old_lambda = as<arma::mat>(out_M["beta_tilde_draw"]);
    old_ssq_m_tilde = as<arma::vec>(out_M["ssq_y_tilde_draw"]);
    //
    oldz_Y = as<arma::vec>(out_Y["zdraw"]);
    oldbeta_2 = as<arma::vec>(out_Y["betadraw"]);
    olddstar_Y = as<arma::mat>(out_Y["dstardraw"]);
    old_y_tilde = as<arma::mat>(out_Y["y_tilde_draw"]);
    old_tau = as<arma::mat>(out_Y["beta_tilde_draw"]);
    old_ssq_y_tilde = as<arma::vec>(out_Y["ssq_y_tilde_draw"]);
    /////////////////////////////////////////////

    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      cutoff_M.slice(mkeep-1) = as<arma::mat>(out_M["cutdraw"]);
      dstardraw_M.slice(mkeep-1) = as<arma::mat>(out_M["dstardraw"]);
      betadraw(mkeep-1,span::all) = trans(as<arma::vec>(out_M["betadraw"]));
      lambdadraw.slice(mkeep-1) = as<arma::mat>(out_M["beta_tilde_draw"]);
      ssq_m_tilde_draw(mkeep-1,span::all) = trans(as<arma::vec>(out_M["ssq_y_tilde_draw"]));
      Mdraw(mkeep-1,span::all) = trans(as<arma::vec>(out_M["zdraw"]));

      cutoff_Y.slice(mkeep-1) = as<arma::mat>(out_Y["cutdraw"]);
      dstardraw_Y.slice(mkeep-1) = as<arma::mat>(out_Y["dstardraw"]);
      beta_2_draw(mkeep-1,span::all) = trans(as<arma::vec>(out_Y["betadraw"]));
      taudraw.slice(mkeep-1) = as<arma::mat>(out_Y["beta_tilde_draw"]);
      ssq_y_tilde_draw(mkeep-1,span::all) = trans(as<arma::vec>(out_Y["ssq_y_tilde_draw"]));
      Ydraw(mkeep-1,span::all) = trans(as<arma::vec>(out_Y["zdraw"]));

      vec mubeta = as<arma::vec>(out_Y["mubeta"]);
      mubeta_2_draw(mkeep-1) = mubeta(nvar_Y-1);  //trans(as<arma::vec>(out_Y["mubeta"]));
      mat varbeta = as<arma::mat>(out_Y["varbeta"]);
      varbeta_2_draw(mkeep-1) = varbeta(nvar_Y-1, nvar_Y-1); //as<arma::mat>(out_Y["varbeta"]);
    }

  }

  return List::create(
    Named("cutoff_M") = cutoff_M,
    // Named("dstardraw_M") = dstardraw_M,
    Named("beta_1") = betadraw,
    Named("lambdadraw") = lambdadraw,
    Named("ssq_m_star_draw") = ssq_m_tilde_draw,
    Named("Mdraw") = Mdraw,
    Named("cutoff_Y") = cutoff_Y,
    // Named("dstardraw_Y") = dstardraw_Y,
    Named("beta_2") = beta_2_draw,
    Named("taudraw") = taudraw,
    Named("ssq_y_star_draw") = ssq_y_tilde_draw,
    Named("Ydraw") = Ydraw,
    Named("mu_draw") = mubeta_2_draw,
    Named("var_draw") = varbeta_2_draw
  );


}
