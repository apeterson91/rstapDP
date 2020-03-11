// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// stapDP_backfit
Rcpp::List stapDP_backfit(const Eigen::VectorXd& y, const Eigen::MatrixXd& Z, const Eigen::VectorXd& X, const double& tau_0, const double& alpha_a, const double& alpha_b, const int& K, const int& iter_max, const int& burn_in, const int& thin, const int& seed, const int& num_posterior_samples);
RcppExport SEXP _rstapDP_stapDP_backfit(SEXP ySEXP, SEXP ZSEXP, SEXP XSEXP, SEXP tau_0SEXP, SEXP alpha_aSEXP, SEXP alpha_bSEXP, SEXP KSEXP, SEXP iter_maxSEXP, SEXP burn_inSEXP, SEXP thinSEXP, SEXP seedSEXP, SEXP num_posterior_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau_0(tau_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_a(alpha_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_b(alpha_bSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_posterior_samples(num_posterior_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(stapDP_backfit(y, Z, X, tau_0, alpha_a, alpha_b, K, iter_max, burn_in, thin, seed, num_posterior_samples));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rstapDP_stapDP_backfit", (DL_FUNC) &_rstapDP_stapDP_backfit, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_rstapDP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
