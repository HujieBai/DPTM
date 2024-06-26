// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// inverse_cpp
Eigen::MatrixXd inverse_cpp(Eigen::MatrixXd& objs);
RcppExport SEXP _DPTM_inverse_cpp(SEXP objsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type objs(objsSEXP);
    rcpp_result_gen = Rcpp::wrap(inverse_cpp(objs));
    return rcpp_result_gen;
END_RCPP
}
// lag_transform
Eigen::MatrixXd lag_transform(Eigen::MatrixXd objs, int t, int n, int lag, bool top);
RcppExport SEXP _DPTM_lag_transform(SEXP objsSEXP, SEXP tSEXP, SEXP nSEXP, SEXP lagSEXP, SEXP topSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type objs(objsSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< bool >::type top(topSEXP);
    rcpp_result_gen = Rcpp::wrap(lag_transform(objs, t, n, lag, top));
    return rcpp_result_gen;
END_RCPP
}
// three_two
double three_two(Eigen::VectorXd& pars1, Eigen::MatrixXd& delty0, Eigen::MatrixXd& evs, Eigen::MatrixXd& omega, int& cd, int& tt, int& nn);
RcppExport SEXP _DPTM_three_two(SEXP pars1SEXP, SEXP delty0SEXP, SEXP evsSEXP, SEXP omegaSEXP, SEXP cdSEXP, SEXP ttSEXP, SEXP nnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type pars1(pars1SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type delty0(delty0SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type evs(evsSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< int& >::type cd(cdSEXP);
    Rcpp::traits::input_parameter< int& >::type tt(ttSEXP);
    Rcpp::traits::input_parameter< int& >::type nn(nnSEXP);
    rcpp_result_gen = Rcpp::wrap(three_two(pars1, delty0, evs, omega, cd, tt, nn));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DPTM_inverse_cpp", (DL_FUNC) &_DPTM_inverse_cpp, 1},
    {"_DPTM_lag_transform", (DL_FUNC) &_DPTM_lag_transform, 5},
    {"_DPTM_three_two", (DL_FUNC) &_DPTM_three_two, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_DPTM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
