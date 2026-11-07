#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <string.h>

/* Kernels for the nodal attribute list of a netsim run. The R callers pass
   only plain integer, double, logical or character vectors without
   attributes (checking attributes needs ATTRIB(), which is not part of R's
   API) and use `[` or c() for everything else, so the result is always what
   base R would have produced. A vector of another type comes back as NULL. */

/* Remove the 1-based positions `drop` from every vector of the list `x`,
   copying the runs between dropped positions with memcpy. Positions outside
   1..length(v) are ignored, as they are by `v[-drop]`; `drop` must hold no NA
   and no zero. Elements of another type come back as NULL. */
SEXP compact_attr_list(SEXP x, SEXP drop) {
  if (TYPEOF(x) != VECSXP) {
    Rf_error("x must be a list");
  }
  PROTECT(drop = Rf_coerceVector(drop, INTSXP));
  R_xlen_t nd_all = XLENGTH(drop);
  int *d = (int *) R_alloc(nd_all, sizeof(int));
  memcpy(d, INTEGER(drop), nd_all * sizeof(int));
  R_isort(d, (int) nd_all);
  for (R_xlen_t k = 0; k < nd_all; k++) {
    if (d[k] == NA_INTEGER || d[k] < 1) {
      Rf_error("drop must hold positive positions without NA");
    }
  }

  R_xlen_t nx = XLENGTH(x);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, nx));
  for (R_xlen_t k = 0; k < nx; k++) {
    SEXP v = VECTOR_ELT(x, k);
    int t = TYPEOF(v);
    if (t != INTSXP && t != REALSXP && t != LGLSXP && t != STRSXP) {
      SET_VECTOR_ELT(out, k, R_NilValue);
      continue;
    }
    R_xlen_t n = XLENGTH(v);
    /* count the distinct positions inside the vector */
    R_xlen_t nd = 0;
    for (R_xlen_t j = 0; j < nd_all; j++) {
      if (d[j] <= n && (j == 0 || d[j] != d[j - 1])) nd++;
    }
    SEXP nv = PROTECT(Rf_allocVector(t, n - nd));
    if (t == STRSXP) {
      R_xlen_t to = 0;
      R_xlen_t j = 0;
      for (R_xlen_t i = 0; i < n; i++) {
        while (j < nd_all && d[j] < i + 1) j++;
        if (j < nd_all && d[j] == i + 1) continue;
        SET_STRING_ELT(nv, to++, STRING_ELT(v, i));
      }
    } else {
      size_t es = (t == REALSXP) ? sizeof(double) : sizeof(int);
      const char *src = (t == REALSXP) ? (const char *) REAL(v)
                        : (t == INTSXP ? (const char *) INTEGER(v)
                                       : (const char *) LOGICAL(v));
      char *dst = (t == REALSXP) ? (char *) REAL(nv)
                  : (t == INTSXP ? (char *) INTEGER(nv) : (char *) LOGICAL(nv));
      R_xlen_t from = 0;
      for (R_xlen_t j = 0; j < nd_all; j++) {
        if (d[j] > n) break;
        R_xlen_t to = d[j] - 1;          /* 0-based position to drop */
        if (to > from) {
          memcpy(dst, src + from * es, (to - from) * es);
          dst += (to - from) * es;
        }
        if (to + 1 > from) from = to + 1; /* duplicates leave `from` as is */
      }
      if (n > from) memcpy(dst, src + from * es, (n - from) * es);
    }
    SET_VECTOR_ELT(out, k, nv);
    UNPROTECT(1);
  }
  Rf_setAttrib(out, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
  UNPROTECT(2);
  return out;
}

/* c(old, new) for attribute-free vectors in the type pairs where c() keeps
   the type of `old`: the same type (integer, double, logical, character),
   integer or logical into double, and logical into integer. NULL otherwise. */
SEXP append_attr_c(SEXP old, SEXP nw) {
  int t = TYPEOF(old), tn = TYPEOF(nw);
  int ok = (t == tn && (t == INTSXP || t == REALSXP || t == LGLSXP || t == STRSXP)) ||
           (t == REALSXP && (tn == INTSXP || tn == LGLSXP)) ||
           (t == INTSXP && tn == LGLSXP);
  if (!ok) {
    return R_NilValue;
  }
  R_xlen_t n1 = XLENGTH(old), n2 = XLENGTH(nw);
  SEXP out = PROTECT(Rf_allocVector(t, n1 + n2));
  if (t == STRSXP) {
    for (R_xlen_t i = 0; i < n1; i++) SET_STRING_ELT(out, i, STRING_ELT(old, i));
    for (R_xlen_t i = 0; i < n2; i++) SET_STRING_ELT(out, n1 + i, STRING_ELT(nw, i));
  } else if (t == REALSXP) {
    double *o = REAL(out);
    memcpy(o, REAL(old), n1 * sizeof(double));
    if (tn == REALSXP) {
      memcpy(o + n1, REAL(nw), n2 * sizeof(double));
    } else {
      const int *b = (tn == INTSXP) ? INTEGER(nw) : LOGICAL(nw);
      for (R_xlen_t i = 0; i < n2; i++) {
        o[n1 + i] = (b[i] == NA_INTEGER) ? NA_REAL : (double) b[i];
      }
    }
  } else {
    int *o = (t == INTSXP) ? INTEGER(out) : LOGICAL(out);
    const int *a = (t == INTSXP) ? INTEGER(old) : LOGICAL(old);
    const int *b = (tn == INTSXP) ? INTEGER(nw) : LOGICAL(nw);
    memcpy(o, a, n1 * sizeof(int));
    memcpy(o + n1, b, n2 * sizeof(int));
  }
  UNPROTECT(1);
  return out;
}
