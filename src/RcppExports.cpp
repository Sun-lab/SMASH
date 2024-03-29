// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Rcpp_ITH_opt
Rcpp::List Rcpp_ITH_opt(const arma::mat& RD, const arma::vec& log_DP, const arma::vec& LBC, const arma::mat& BB, const arma::umat& uniq_BB, const arma::mat& uniq_CN_MA, const arma::mat& eS, const double& purity, const arma::vec& tCN, const double& pi_eps0, const arma::vec& pi0, const arma::vec& unc_qq0, const bool& mstep, const arma::uword& max_iter, const double& eps, const bool& show);
RcppExport SEXP _SMASH_Rcpp_ITH_opt(SEXP RDSEXP, SEXP log_DPSEXP, SEXP LBCSEXP, SEXP BBSEXP, SEXP uniq_BBSEXP, SEXP uniq_CN_MASEXP, SEXP eSSEXP, SEXP puritySEXP, SEXP tCNSEXP, SEXP pi_eps0SEXP, SEXP pi0SEXP, SEXP unc_qq0SEXP, SEXP mstepSEXP, SEXP max_iterSEXP, SEXP epsSEXP, SEXP showSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type RD(RDSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type log_DP(log_DPSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type LBC(LBCSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type BB(BBSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type uniq_BB(uniq_BBSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type uniq_CN_MA(uniq_CN_MASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type eS(eSSEXP);
    Rcpp::traits::input_parameter< const double& >::type purity(puritySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tCN(tCNSEXP);
    Rcpp::traits::input_parameter< const double& >::type pi_eps0(pi_eps0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi0(pi0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type unc_qq0(unc_qq0SEXP);
    Rcpp::traits::input_parameter< const bool& >::type mstep(mstepSEXP);
    Rcpp::traits::input_parameter< const arma::uword& >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show(showSEXP);
    rcpp_result_gen = Rcpp::wrap(Rcpp_ITH_opt(RD, log_DP, LBC, BB, uniq_BB, uniq_CN_MA, eS, purity, tCN, pi_eps0, pi0, unc_qq0, mstep, max_iter, eps, show));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SMASH_Rcpp_ITH_opt", (DL_FUNC) &_SMASH_Rcpp_ITH_opt, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_SMASH(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
