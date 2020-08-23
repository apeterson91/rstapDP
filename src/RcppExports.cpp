// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// stappDP_fit
Rcpp::List stappDP_fit(const Eigen::VectorXd& y, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X, const Eigen::MatrixXd& S, const Eigen::VectorXd& w, const double& alpha_a, const double& alpha_b, const double& sigma_a, const double& sigma_b, const double& tau_a, const double& tau_b, const int& K, const int& num_penalties, const int& iter_max, const int& burn_in, const int& thin, const int& seed, const int& num_posterior_samples, const bool& fix_alpha);
RcppExport SEXP _rstapDP_stappDP_fit(SEXP ySEXP, SEXP ZSEXP, SEXP XSEXP, SEXP SSEXP, SEXP wSEXP, SEXP alpha_aSEXP, SEXP alpha_bSEXP, SEXP sigma_aSEXP, SEXP sigma_bSEXP, SEXP tau_aSEXP, SEXP tau_bSEXP, SEXP KSEXP, SEXP num_penaltiesSEXP, SEXP iter_maxSEXP, SEXP burn_inSEXP, SEXP thinSEXP, SEXP seedSEXP, SEXP num_posterior_samplesSEXP, SEXP fix_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_a(alpha_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_b(alpha_bSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma_a(sigma_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma_b(sigma_bSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_a(tau_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_b(tau_bSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_penalties(num_penaltiesSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_posterior_samples(num_posterior_samplesSEXP);
    Rcpp::traits::input_parameter< const bool& >::type fix_alpha(fix_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(stappDP_fit(y, Z, X, S, w, alpha_a, alpha_b, sigma_a, sigma_b, tau_a, tau_b, K, num_penalties, iter_max, burn_in, thin, seed, num_posterior_samples, fix_alpha));
    return rcpp_result_gen;
END_RCPP
}
// stappDP_mer_fit
Rcpp::List stappDP_mer_fit(const Eigen::VectorXd& y, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X, const Eigen::ArrayXXd& W, const Eigen::MatrixXd& S, const Eigen::VectorXd& w, const SEXP& subj_mat_, const Eigen::ArrayXi& subj_n, const double& alpha_a, const double& alpha_b, const double& sigma_a, const double& sigma_b, const double& tau_a, const double& tau_b, const int& K, const int& num_penalties, const int& iter_max, const int& burn_in, const int& thin, const int& seed, const int& num_posterior_samples, const bool& fix_alpha);
RcppExport SEXP _rstapDP_stappDP_mer_fit(SEXP ySEXP, SEXP ZSEXP, SEXP XSEXP, SEXP WSEXP, SEXP SSEXP, SEXP wSEXP, SEXP subj_mat_SEXP, SEXP subj_nSEXP, SEXP alpha_aSEXP, SEXP alpha_bSEXP, SEXP sigma_aSEXP, SEXP sigma_bSEXP, SEXP tau_aSEXP, SEXP tau_bSEXP, SEXP KSEXP, SEXP num_penaltiesSEXP, SEXP iter_maxSEXP, SEXP burn_inSEXP, SEXP thinSEXP, SEXP seedSEXP, SEXP num_posterior_samplesSEXP, SEXP fix_alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const SEXP& >::type subj_mat_(subj_mat_SEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXi& >::type subj_n(subj_nSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_a(alpha_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_b(alpha_bSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma_a(sigma_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma_b(sigma_bSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_a(tau_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_b(tau_bSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_penalties(num_penaltiesSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_posterior_samples(num_posterior_samplesSEXP);
    Rcpp::traits::input_parameter< const bool& >::type fix_alpha(fix_alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(stappDP_mer_fit(y, Z, X, W, S, w, subj_mat_, subj_n, alpha_a, alpha_b, sigma_a, sigma_b, tau_a, tau_b, K, num_penalties, iter_max, burn_in, thin, seed, num_posterior_samples, fix_alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rstapDP_stappDP_fit", (DL_FUNC) &_rstapDP_stappDP_fit, 19},
    {"_rstapDP_stappDP_mer_fit", (DL_FUNC) &_rstapDP_stappDP_mer_fit, 22},
    {NULL, NULL, 0}
};

RcppExport void R_init_rstapDP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
