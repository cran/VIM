// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gowerd
RcppExport SEXP gowerd(SEXP dataX, SEXP dataY, SEXP weights, SEXP ncolNUMFAC, SEXP levOrders, SEXP mixedConstants);
RcppExport SEXP _VIM_gowerd(SEXP dataXSEXP, SEXP dataYSEXP, SEXP weightsSEXP, SEXP ncolNUMFACSEXP, SEXP levOrdersSEXP, SEXP mixedConstantsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type dataX(dataXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dataY(dataYSEXP);
    Rcpp::traits::input_parameter< SEXP >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ncolNUMFAC(ncolNUMFACSEXP);
    Rcpp::traits::input_parameter< SEXP >::type levOrders(levOrdersSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mixedConstants(mixedConstantsSEXP);
    rcpp_result_gen = Rcpp::wrap(gowerd(dataX, dataY, weights, ncolNUMFAC, levOrders, mixedConstants));
    return rcpp_result_gen;
END_RCPP
}
// whichminN
RcppExport SEXP whichminN(SEXP xR, SEXP nR, int returnValue);
RcppExport SEXP _VIM_whichminN(SEXP xRSEXP, SEXP nRSEXP, SEXP returnValueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xR(xRSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nR(nRSEXP);
    Rcpp::traits::input_parameter< int >::type returnValue(returnValueSEXP);
    rcpp_result_gen = Rcpp::wrap(whichminN(xR, nR, returnValue));
    return rcpp_result_gen;
END_RCPP
}
// gowerDind
RcppExport SEXP gowerDind(SEXP dataX, SEXP dataY, SEXP weights, SEXP ncolNUMFAC, SEXP levOrders, SEXP mixedConstants, SEXP nR, SEXP returnMinR);
RcppExport SEXP _VIM_gowerDind(SEXP dataXSEXP, SEXP dataYSEXP, SEXP weightsSEXP, SEXP ncolNUMFACSEXP, SEXP levOrdersSEXP, SEXP mixedConstantsSEXP, SEXP nRSEXP, SEXP returnMinRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type dataX(dataXSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dataY(dataYSEXP);
    Rcpp::traits::input_parameter< SEXP >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ncolNUMFAC(ncolNUMFACSEXP);
    Rcpp::traits::input_parameter< SEXP >::type levOrders(levOrdersSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mixedConstants(mixedConstantsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nR(nRSEXP);
    Rcpp::traits::input_parameter< SEXP >::type returnMinR(returnMinRSEXP);
    rcpp_result_gen = Rcpp::wrap(gowerDind(dataX, dataY, weights, ncolNUMFAC, levOrders, mixedConstants, nR, returnMinR));
    return rcpp_result_gen;
END_RCPP
}
