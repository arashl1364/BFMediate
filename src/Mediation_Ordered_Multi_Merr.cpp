#include "BFMediate.h"

// // [[Rcpp::depends(RcppArmadillo)]]
// // #include "bayesm.h" //if you include bayesm.h here, you do not need to define functions before calling them
// #include "RcppArmadillo.h"
// // #include "rordprobitGibbs_me.cpp"
// // #include "rordprobitGibbs_me_M.cpp"
// // #include "rordprobitGibbs_me_MY.cpp"
// using namespace arma; // use the Armadillo library for matrix computations
// using namespace Rcpp;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////            MAIN FUNCTION                  /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
  arma::vec M = randu<arma::vec>(ny);
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
  arma::vec old_ssq_y_tilde(Y_ind);
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
