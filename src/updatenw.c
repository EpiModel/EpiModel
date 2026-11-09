#include <R.h>
#include <Rinternals.h>

/* Shift sorted vertex ids down past deleted ids.

   orig: vertex ids sorted in increasing order
   del:  deleted vertex ids sorted in increasing order

   Returns orig[pos_orig] minus the number of elements of del that are
   <= orig[pos_orig], as an integer vector, so that the ids left after
   deleting del are consecutive again. Both inputs are coerced to integer.

   Example: orig = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 and del = 4, 6, 8, 9 return
   1, 2, 3, 3, 4, 4, 5, 5, 5, 6, so the ids kept (1, 2, 3, 5, 7, 10) become
   1 to 6. */
SEXP shiftVec(SEXP orig, SEXP del) {
  if (!Rf_isNumeric(orig) || !Rf_isNumeric(del)) {
    Rf_error("orig and del must be numeric vectors");
  }
  PROTECT(orig = Rf_coerceVector(orig, INTSXP));
  PROTECT(del = Rf_coerceVector(del, INTSXP));
  R_xlen_t n_orig = XLENGTH(orig);
  R_xlen_t n_del = XLENGTH(del);
  int *p_orig = INTEGER(orig);
  int *p_del = INTEGER(del);

  SEXP shifted = PROTECT(Rf_allocVector(INTSXP, n_orig));
  int *p_shifted = INTEGER(shifted);

  /* the ids in orig from del[pos_del - 1] up to but not including
     del[pos_del] have exactly pos_del deleted ids at or below them */
  R_xlen_t pos_orig = 0;
  for (R_xlen_t pos_del = 0; pos_del < n_del; pos_del++) {
    while (pos_orig < n_orig && p_orig[pos_orig] < p_del[pos_del]) {
      p_shifted[pos_orig] = p_orig[pos_orig] - (int) pos_del;
      pos_orig++;
    }
  }

  /* any trailing part of orig after the end of del */
  while (pos_orig < n_orig) {
    p_shifted[pos_orig] = p_orig[pos_orig] - (int) n_del;
    pos_orig++;
  }

  UNPROTECT(3);
  return shifted;
}
