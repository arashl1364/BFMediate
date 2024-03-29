// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// runiregGibbs_rcpp_me
List runiregGibbs_rcpp_me(arma::vec const& y, arma::mat const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq, double sigmasq, int R, int keep, int nprint, bool betafix, bool sigmafix, arma::mat betavalue, arma::vec sigmavalue);
RcppExport SEXP _BFMediate_runiregGibbs_rcpp_me(SEXP ySEXP, SEXP XSEXP, SEXP betabarSEXP, SEXP ASEXP, SEXP nuSEXP, SEXP ssqSEXP, SEXP sigmasqSEXP, SEXP RSEXP, SEXP keepSEXP, SEXP nprintSEXP, SEXP betafixSEXP, SEXP sigmafixSEXP, SEXP betavalueSEXP, SEXP sigmavalueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betabar(betabarSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type ssq(ssqSEXP);
    Rcpp::traits::input_parameter< double >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< int >::type nprint(nprintSEXP);
    Rcpp::traits::input_parameter< bool >::type betafix(betafixSEXP);
    Rcpp::traits::input_parameter< bool >::type sigmafix(sigmafixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type betavalue(betavalueSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigmavalue(sigmavalueSEXP);
    rcpp_result_gen = Rcpp::wrap(runiregGibbs_rcpp_me(y, X, betabar, A, nu, ssq, sigmasq, R, keep, nprint, betafix, sigmafix, betavalue, sigmavalue));
    return rcpp_result_gen;
END_RCPP
}
// MeasurementMCatCpp
List MeasurementMCatCpp(arma::vec const& dep, arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad, arma::mat const& A_2, arma::vec const& betabar_2, double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat, int const& Y_ind, int R, int keep, int nprint);
RcppExport SEXP _BFMediate_MeasurementMCatCpp(SEXP depSEXP, SEXP ySEXP, SEXP XSEXP, SEXP kSEXP, SEXP ASEXP, SEXP betabarSEXP, SEXP AdSEXP, SEXP A_2SEXP, SEXP betabar_2SEXP, SEXP sSEXP, SEXP inc_rootSEXP, SEXP dstarbarSEXP, SEXP betahatSEXP, SEXP Y_indSEXP, SEXP RSEXP, SEXP keepSEXP, SEXP nprintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec const& >::type dep(depSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betabar(betabarSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Ad(AdSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type A_2(A_2SEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betabar_2(betabar_2SEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type inc_root(inc_rootSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type dstarbar(dstarbarSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< int const& >::type Y_ind(Y_indSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< int >::type nprint(nprintSEXP);
    rcpp_result_gen = Rcpp::wrap(MeasurementMCatCpp(dep, y, X, k, A, betabar, Ad, A_2, betabar_2, s, inc_root, dstarbar, betahat, Y_ind, R, keep, nprint));
    return rcpp_result_gen;
END_RCPP
}
// MeasurementYCatCpp
List MeasurementYCatCpp(arma::mat const& y, arma::mat const& X, int k, arma::mat const& A, arma::vec const& betabar, arma::mat const& Ad, double s, arma::mat const& inc_root, arma::vec const& dstarbar, arma::vec const& betahat, int const& Y_ind, int R, int keep, int nprint);
RcppExport SEXP _BFMediate_MeasurementYCatCpp(SEXP ySEXP, SEXP XSEXP, SEXP kSEXP, SEXP ASEXP, SEXP betabarSEXP, SEXP AdSEXP, SEXP sSEXP, SEXP inc_rootSEXP, SEXP dstarbarSEXP, SEXP betahatSEXP, SEXP Y_indSEXP, SEXP RSEXP, SEXP keepSEXP, SEXP nprintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betabar(betabarSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Ad(AdSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type inc_root(inc_rootSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type dstarbar(dstarbarSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< int const& >::type Y_ind(Y_indSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< int >::type nprint(nprintSEXP);
    rcpp_result_gen = Rcpp::wrap(MeasurementYCatCpp(y, X, k, A, betabar, Ad, s, inc_root, dstarbar, betahat, Y_ind, R, keep, nprint));
    return rcpp_result_gen;
END_RCPP
}
// MeasurementMYCatCpp
List MeasurementMYCatCpp(arma::mat const& X, arma::mat const& m_star, arma::mat const& y_star, int k_M, int k_Y, int M_ind, int Y_ind, arma::mat const& A_M, arma::vec const& betabar, arma::mat const& Ad_M, double s_M, arma::mat const& inc_root_M, arma::vec const& dstarbar_M, arma::vec const& betahat, arma::mat const& A_Y, arma::vec const& beta_2_bar, arma::mat const& Ad_Y, double s_Y, arma::mat const& inc_root_Y, arma::vec const& dstarbar_Y, arma::vec const& beta_2_hat, int R, int keep, int nprint);
RcppExport SEXP _BFMediate_MeasurementMYCatCpp(SEXP XSEXP, SEXP m_starSEXP, SEXP y_starSEXP, SEXP k_MSEXP, SEXP k_YSEXP, SEXP M_indSEXP, SEXP Y_indSEXP, SEXP A_MSEXP, SEXP betabarSEXP, SEXP Ad_MSEXP, SEXP s_MSEXP, SEXP inc_root_MSEXP, SEXP dstarbar_MSEXP, SEXP betahatSEXP, SEXP A_YSEXP, SEXP beta_2_barSEXP, SEXP Ad_YSEXP, SEXP s_YSEXP, SEXP inc_root_YSEXP, SEXP dstarbar_YSEXP, SEXP beta_2_hatSEXP, SEXP RSEXP, SEXP keepSEXP, SEXP nprintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type m_star(m_starSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type y_star(y_starSEXP);
    Rcpp::traits::input_parameter< int >::type k_M(k_MSEXP);
    Rcpp::traits::input_parameter< int >::type k_Y(k_YSEXP);
    Rcpp::traits::input_parameter< int >::type M_ind(M_indSEXP);
    Rcpp::traits::input_parameter< int >::type Y_ind(Y_indSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type A_M(A_MSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betabar(betabarSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Ad_M(Ad_MSEXP);
    Rcpp::traits::input_parameter< double >::type s_M(s_MSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type inc_root_M(inc_root_MSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type dstarbar_M(dstarbar_MSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type A_Y(A_YSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type beta_2_bar(beta_2_barSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Ad_Y(Ad_YSEXP);
    Rcpp::traits::input_parameter< double >::type s_Y(s_YSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type inc_root_Y(inc_root_YSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type dstarbar_Y(dstarbar_YSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type beta_2_hat(beta_2_hatSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< int >::type nprint(nprintSEXP);
    rcpp_result_gen = Rcpp::wrap(MeasurementMYCatCpp(X, m_star, y_star, k_M, k_Y, M_ind, Y_ind, A_M, betabar, Ad_M, s_M, inc_root_M, dstarbar_M, betahat, A_Y, beta_2_bar, Ad_Y, s_Y, inc_root_Y, dstarbar_Y, beta_2_hat, R, keep, nprint));
    return rcpp_result_gen;
END_RCPP
}
// RuniregGibbsMultiCpp
List RuniregGibbsMultiCpp(arma::mat const& y, arma::mat const& M, arma::vec const& X, arma::vec const& betabar, arma::mat const& A, double nu, double ssq, arma::vec sigmasq, int R, int keep, int nprint, bool betafix, bool sigmafix, arma::mat betavalue, arma::vec sigmavalue);
RcppExport SEXP _BFMediate_RuniregGibbsMultiCpp(SEXP ySEXP, SEXP MSEXP, SEXP XSEXP, SEXP betabarSEXP, SEXP ASEXP, SEXP nuSEXP, SEXP ssqSEXP, SEXP sigmasqSEXP, SEXP RSEXP, SEXP keepSEXP, SEXP nprintSEXP, SEXP betafixSEXP, SEXP sigmafixSEXP, SEXP betavalueSEXP, SEXP sigmavalueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type betabar(betabarSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type ssq(ssqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< int >::type nprint(nprintSEXP);
    Rcpp::traits::input_parameter< bool >::type betafix(betafixSEXP);
    Rcpp::traits::input_parameter< bool >::type sigmafix(sigmafixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type betavalue(betavalueSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigmavalue(sigmavalueSEXP);
    rcpp_result_gen = Rcpp::wrap(RuniregGibbsMultiCpp(y, M, X, betabar, A, nu, ssq, sigmasq, R, keep, nprint, betafix, sigmafix, betavalue, sigmavalue));
    return rcpp_result_gen;
END_RCPP
}
// rbprobitGibbs_rcpp_me
List rbprobitGibbs_rcpp_me(arma::vec const& y, arma::mat const& X, arma::vec const& Abetabar, arma::mat const& root, arma::vec beta, arma::vec const& sigma, arma::vec const& trunpt, arma::vec const& above, int R, int keep, int nprint);
RcppExport SEXP _BFMediate_rbprobitGibbs_rcpp_me(SEXP ySEXP, SEXP XSEXP, SEXP AbetabarSEXP, SEXP rootSEXP, SEXP betaSEXP, SEXP sigmaSEXP, SEXP trunptSEXP, SEXP aboveSEXP, SEXP RSEXP, SEXP keepSEXP, SEXP nprintSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type Abetabar(AbetabarSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type root(rootSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type trunpt(trunptSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type above(aboveSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type keep(keepSEXP);
    Rcpp::traits::input_parameter< int >::type nprint(nprintSEXP);
    rcpp_result_gen = Rcpp::wrap(rbprobitGibbs_rcpp_me(y, X, Abetabar, root, beta, sigma, trunpt, above, R, keep, nprint));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BFMediate_runiregGibbs_rcpp_me", (DL_FUNC) &_BFMediate_runiregGibbs_rcpp_me, 14},
    {"_BFMediate_MeasurementMCatCpp", (DL_FUNC) &_BFMediate_MeasurementMCatCpp, 17},
    {"_BFMediate_MeasurementYCatCpp", (DL_FUNC) &_BFMediate_MeasurementYCatCpp, 14},
    {"_BFMediate_MeasurementMYCatCpp", (DL_FUNC) &_BFMediate_MeasurementMYCatCpp, 24},
    {"_BFMediate_RuniregGibbsMultiCpp", (DL_FUNC) &_BFMediate_RuniregGibbsMultiCpp, 15},
    {"_BFMediate_rbprobitGibbs_rcpp_me", (DL_FUNC) &_BFMediate_rbprobitGibbs_rcpp_me, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_BFMediate(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
