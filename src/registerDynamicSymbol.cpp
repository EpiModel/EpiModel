// RegisteringDynamic Symbols

#include <Rcpp.h>
using namespace Rcpp;

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// [[Rcpp::init]]
void EpiModel_init(DllInfo* info) {
  R_useDynamicSymbols(info, TRUE);
}
