// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dap
List dap(List arg);
RcppExport SEXP _dap_dap(SEXP argSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type arg(argSEXP);
    rcpp_result_gen = Rcpp::wrap(dap(arg));
    return rcpp_result_gen;
END_RCPP
}
// read_sbams
List read_sbams(const char* file, bool normalize);
RcppExport SEXP _dap_read_sbams(SEXP fileSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type file(fileSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(read_sbams(file, normalize));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dap_dap", (DL_FUNC) &_dap_dap, 1},
    {"_dap_read_sbams", (DL_FUNC) &_dap_read_sbams, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_dap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
