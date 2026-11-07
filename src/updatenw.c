#include <R.h>
#include <Rinternals.h>

/* Shift sorted vertex ids down past deleted ids.

   x: vertex ids sorted in increasing order
   y: deleted vertex ids sorted in increasing order

   Returns x[i] minus the number of elements of y that are <= x[i], as an
   integer vector, so that the ids left after deleting y are consecutive
   again. Both inputs are coerced to integer. */
SEXP shiftVec(SEXP x, SEXP y) {
  if (!Rf_isNumeric(x) || !Rf_isNumeric(y)) {
    Rf_error("x and y must be numeric vectors");
  }
  PROTECT(x = Rf_coerceVector(x, INTSXP));
  PROTECT(y = Rf_coerceVector(y, INTSXP));
  R_xlen_t x_len = XLENGTH(x);
  R_xlen_t y_len = XLENGTH(y);
  int *px = INTEGER(x);
  int *py = INTEGER(y);

  SEXP z = PROTECT(Rf_allocVector(INTSXP, x_len));
  int *pz = INTEGER(z);

  R_xlen_t i = 0;
  for (R_xlen_t j = 0; j < y_len; j++) {
    while (i < x_len && px[i] < py[j]) {
      pz[i] = px[i] - (int) j;
      i++;
    }
  }

  /* any trailing part of x after the end of y */
  while (i < x_len) {
    pz[i] = px[i] - (int) y_len;
    i++;
  }

  UNPROTECT(3);
  return z;
}
