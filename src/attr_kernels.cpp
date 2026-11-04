#include <Rcpp.h>
#include <cstring>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Remove the (1-based) positions `drop` from every vector of the attribute
// list `x` in one pass per vector, copying the runs between the dropped
// positions with memcpy. `drop` must hold valid positions (1 to length).
// Elements that are not plain integer, double, or logical vectors (or that
// carry R attributes) come back as NULL, for the R caller to subset.
// [[Rcpp::export(rng = false)]]
List compact_attr_list_cpp(List x, IntegerVector drop) {
  std::vector<int> d(drop.begin(), drop.end());
  std::sort(d.begin(), d.end());
  d.erase(std::unique(d.begin(), d.end()), d.end());
  R_xlen_t nx = x.size();
  List out(nx);
  for (R_xlen_t k = 0; k < nx; k++) {
    SEXP v = x[k];
    int t = TYPEOF(v);
    if (ANY_ATTRIB(v) || (t != INTSXP && t != REALSXP && t != LGLSXP)) {
      out[k] = R_NilValue;
      continue;
    }
    R_xlen_t n = XLENGTH(v);
    size_t es = (t == REALSXP) ? sizeof(double) : sizeof(int);
    R_xlen_t nd = 0;
    for (int di : d) if (di >= 1 && di <= n) nd++;
    SEXP nv = PROTECT(Rf_allocVector(t, n - nd));
    const char* src = (t == REALSXP) ? (const char*) REAL(v) :
      (t == INTSXP ? (const char*) INTEGER(v) : (const char*) LOGICAL(v));
    char* dst = (t == REALSXP) ? (char*) REAL(nv) :
      (t == INTSXP ? (char*) INTEGER(nv) : (char*) LOGICAL(nv));
    R_xlen_t from = 0;
    for (int di : d) {
      if (di < 1 || di > n) continue;
      R_xlen_t to = di - 1;
      if (to > from) {
        std::memcpy(dst, src + from * es, (to - from) * es);
        dst += (to - from) * es;
      }
      from = di;
    }
    if (n > from) std::memcpy(dst, src + from * es, (n - from) * es);
    out[k] = nv;
    UNPROTECT(1);
  }
  out.attr("names") = x.attr("names");
  return out;
}

// c(old, new) by memcpy, for attribute-free atomic vectors and the type pairs
// where c() keeps the type of `old` (same type; integer or logical into
// double; logical into integer). Returns NULL otherwise, and the R caller then
// uses c().
// [[Rcpp::export(rng = false)]]
SEXP append_attr_cpp(SEXP old, SEXP nw) {
  int t = TYPEOF(old), tn = TYPEOF(nw);
  if (ANY_ATTRIB(old) || ANY_ATTRIB(nw)) return R_NilValue;
  bool ok = (t == tn && (t == INTSXP || t == REALSXP || t == LGLSXP)) ||
            (t == REALSXP && (tn == INTSXP || tn == LGLSXP)) ||
            (t == INTSXP && tn == LGLSXP);
  if (!ok) return R_NilValue;
  R_xlen_t n1 = XLENGTH(old), n2 = XLENGTH(nw);
  SEXP out = PROTECT(Rf_allocVector(t, n1 + n2));
  if (t == REALSXP) {
    double* o = REAL(out);
    std::memcpy(o, REAL(old), n1 * sizeof(double));
    if (tn == REALSXP) {
      std::memcpy(o + n1, REAL(nw), n2 * sizeof(double));
    } else {
      const int* b = (tn == INTSXP) ? INTEGER(nw) : LOGICAL(nw);
      for (R_xlen_t i = 0; i < n2; i++)
        o[n1 + i] = (b[i] == NA_INTEGER) ? NA_REAL : (double) b[i];
    }
  } else {
    int* o = (t == INTSXP) ? INTEGER(out) : LOGICAL(out);
    const int* a = (t == INTSXP) ? INTEGER(old) : LOGICAL(old);
    const int* b = (tn == INTSXP) ? INTEGER(nw) : LOGICAL(nw);
    std::memcpy(o, a, n1 * sizeof(int));
    std::memcpy(o + n1, b, n2 * sizeof(int));
  }
  UNPROTECT(1);
  return out;
}
