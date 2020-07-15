// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// stappDP_fit
Rcpp::List stappDP_fit(const Eigen::VectorXd& y, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X, const Eigen::MatrixXd& S, const Eigen::VectorXd& w, const double& nu_0, const double& alpha_a, const double& alpha_b, const int& K, const int& num_penalties, const int& iter_max, const int& burn_in, const int& thin, const int& seed, const int& num_posterior_samples);
RcppExport SEXP _rstapDP_stappDP_fit(SEXP ySEXP, SEXP ZSEXP, SEXP XSEXP, SEXP SSEXP, SEXP wSEXP, SEXP nu_0SEXP, SEXP alpha_aSEXP, SEXP alpha_bSEXP, SEXP KSEXP, SEXP num_penaltiesSEXP, SEXP iter_maxSEXP, SEXP burn_inSEXP, SEXP thinSEXP, SEXP seedSEXP, SEXP num_posterior_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double& >::type nu_0(nu_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_a(alpha_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_b(alpha_bSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_penalties(num_penaltiesSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_posterior_samples(num_posterior_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(stappDP_fit(y, Z, X, S, w, nu_0, alpha_a, alpha_b, K, num_penalties, iter_max, burn_in, thin, seed, num_posterior_samples));
    return rcpp_result_gen;
END_RCPP
}
// stapDP_fit
Rcpp::List stapDP_fit(const Eigen::VectorXd& y, const Eigen::MatrixXd& Z, const Eigen::MatrixXd& X, const Eigen::VectorXd& w, const Eigen::ArrayXd tau_0, const double& alpha_a, const double& alpha_b, const int& K, const int& iter_max, const int& burn_in, const int& thin, const int& seed, const int& num_posterior_samples);
RcppExport SEXP _rstapDP_stapDP_fit(SEXP ySEXP, SEXP ZSEXP, SEXP XSEXP, SEXP wSEXP, SEXP tau_0SEXP, SEXP alpha_aSEXP, SEXP alpha_bSEXP, SEXP KSEXP, SEXP iter_maxSEXP, SEXP burn_inSEXP, SEXP thinSEXP, SEXP seedSEXP, SEXP num_posterior_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXd >::type tau_0(tau_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_a(alpha_aSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_b(alpha_bSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter_max(iter_maxSEXP);
    Rcpp::traits::input_parameter< const int& >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_posterior_samples(num_posterior_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(stapDP_fit(y, Z, X, w, tau_0, alpha_a, alpha_b, K, iter_max, burn_in, thin, seed, num_posterior_samples));
    return rcpp_result_gen;
END_RCPP
}
// green_loss_unknown
Eigen::ArrayXd green_loss_unknown(const Eigen::ArrayXXi& cluster_assignment, const Eigen::ArrayXXd& pmat, const double& tau);
RcppExport SEXP _rstapDP_green_loss_unknown(SEXP cluster_assignmentSEXP, SEXP pmatSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXi& >::type cluster_assignment(cluster_assignmentSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type pmat(pmatSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(green_loss_unknown(cluster_assignment, pmat, tau));
    return rcpp_result_gen;
END_RCPP
}
// green_loss_known
Eigen::ArrayXd green_loss_known(const Eigen::ArrayXXi& cluster_assignment, const Eigen::ArrayXXd& pmat, const Eigen::ArrayXXi& true_cluster_assignment, const double& a, const double& b);
RcppExport SEXP _rstapDP_green_loss_known(SEXP cluster_assignmentSEXP, SEXP pmatSEXP, SEXP true_cluster_assignmentSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXi& >::type cluster_assignment(cluster_assignmentSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type pmat(pmatSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXi& >::type true_cluster_assignment(true_cluster_assignmentSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(green_loss_known(cluster_assignment, pmat, true_cluster_assignment, a, b));
    return rcpp_result_gen;
END_RCPP
}
// square_error
Eigen::ArrayXd square_error(const Eigen::ArrayXXi& cluster_assignment, const Eigen::ArrayXXd& pmat);
RcppExport SEXP _rstapDP_square_error(SEXP cluster_assignmentSEXP, SEXP pmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXi& >::type cluster_assignment(cluster_assignmentSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type pmat(pmatSEXP);
    rcpp_result_gen = Rcpp::wrap(square_error(cluster_assignment, pmat));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rstapDP_stappDP_fit", (DL_FUNC) &_rstapDP_stappDP_fit, 15},
    {"_rstapDP_stapDP_fit", (DL_FUNC) &_rstapDP_stapDP_fit, 13},
    {"_rstapDP_green_loss_unknown", (DL_FUNC) &_rstapDP_green_loss_unknown, 3},
    {"_rstapDP_green_loss_known", (DL_FUNC) &_rstapDP_green_loss_known, 5},
    {"_rstapDP_square_error", (DL_FUNC) &_rstapDP_square_error, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rstapDP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
