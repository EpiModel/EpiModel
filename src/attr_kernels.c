#include <R.h>
#include <Rinternals.h>
#include <string.h>

/* Kernels for the nodal attribute list of a netsim run, called from
   remove_node_attr() and append_attr() in R/net.fn.accessor.R. Those callers
   validate the inputs and pass only plain integer, double, logical or
   character vectors without attributes (checking attributes needs ATTRIB(),
   which is not part of R's API); everything else goes through `[` or c() in
   R, so the result is always what base R would have produced. Any other type
   reaching these functions is an error. */

/* Copy the `n_elem` elements of `src`, each `elem_size` bytes, to `dst`,
   skipping the 1-based positions drop_pos[0], ..., drop_pos[n_drop - 1],
   which are increasing and at most `n_elem`. Each run of kept elements
   between two dropped positions is one memcpy. */
static void copy_skipping(char *dst, const char *src, R_xlen_t n_elem,
                          size_t elem_size, const int *drop_pos,
                          R_xlen_t n_drop) {
  R_xlen_t run_start = 0; /* 0-based start of the current run of kept elements */
  for (R_xlen_t k = 0; k < n_drop; k++) {
    R_xlen_t run_end = drop_pos[k] - 1; /* 0-based dropped position */
    R_xlen_t run_len = run_end - run_start;
    if (run_len > 0) {
      memcpy(dst, src + run_start * elem_size, run_len * elem_size);
      dst += run_len * elem_size;
    }
    run_start = run_end + 1;
  }
  if (n_elem > run_start) {
    memcpy(dst, src + run_start * elem_size, (n_elem - run_start) * elem_size);
  }
}

/* Remove the positions `drop_pos` from every vector of the list `attr_list`.

   attr_list: list of plain integer, double, logical or character vectors
   drop_pos:  1-based positions to drop, an integer vector in increasing order
              with no duplicates, NA or values below 1 (remove_node_attr()
              validates the positions and passes sort(unique(.)))

   Positions past the end of a vector are ignored, as they are by
   `v[-drop_pos]`. Returns a list of the compacted vectors, with the names of
   `attr_list`. */
SEXP compact_attr_list(SEXP attr_list, SEXP drop_pos) {
  if (TYPEOF(attr_list) != VECSXP || TYPEOF(drop_pos) != INTSXP) {
    Rf_error("attr_list must be a list and drop_pos an integer vector");
  }
  const int *drop = INTEGER(drop_pos);
  R_xlen_t n_drop = XLENGTH(drop_pos);
  R_xlen_t n_attr = XLENGTH(attr_list);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, n_attr));
  for (R_xlen_t a = 0; a < n_attr; a++) {
    SEXP old_vec = VECTOR_ELT(attr_list, a);
    int type = TYPEOF(old_vec);
    if (type != INTSXP && type != REALSXP && type != LGLSXP && type != STRSXP) {
      Rf_error("element %d of attr_list is not an integer, double, logical "
               "or character vector", (int) (a + 1));
    }
    R_xlen_t old_len = XLENGTH(old_vec);

    /* drop_pos is increasing, so the positions inside this vector come first */
    R_xlen_t n_drop_in = n_drop;
    while (n_drop_in > 0 && drop[n_drop_in - 1] > old_len) n_drop_in--;

    SEXP new_vec = Rf_allocVector(type, old_len - n_drop_in);
    SET_VECTOR_ELT(out, a, new_vec);

    if (type == STRSXP) {
      R_xlen_t next_pos = 0; /* next position of new_vec to fill */
      R_xlen_t k = 0;        /* index in drop of the next position to skip */
      for (R_xlen_t pos = 0; pos < old_len; pos++) {
        if (k < n_drop_in && drop[k] == pos + 1) {
          k++;
          continue;
        }
        SET_STRING_ELT(new_vec, next_pos++, STRING_ELT(old_vec, pos));
      }
    } else if (type == REALSXP) {
      copy_skipping((char *) REAL(new_vec), (const char *) REAL(old_vec),
                    old_len, sizeof(double), drop, n_drop_in);
    } else if (type == INTSXP) {
      copy_skipping((char *) INTEGER(new_vec), (const char *) INTEGER(old_vec),
                    old_len, sizeof(int), drop, n_drop_in);
    } else {
      copy_skipping((char *) LOGICAL(new_vec), (const char *) LOGICAL(old_vec),
                    old_len, sizeof(int), drop, n_drop_in);
    }
  }
  Rf_setAttrib(out, R_NamesSymbol, Rf_getAttrib(attr_list, R_NamesSymbol));
  UNPROTECT(1);
  return out;
}

/* c(old_vec, new_vals) for two plain vectors in the type pairs where c() keeps
   the type of `old_vec`: the same type (integer, double, logical, character),
   integer or logical into double, and logical into integer. append_attr()
   checks the pair before calling; any other pair is an error. */
SEXP append_attr_c(SEXP old_vec, SEXP new_vals) {
  int old_type = TYPEOF(old_vec);
  int new_type = TYPEOF(new_vals);
  R_xlen_t old_len = XLENGTH(old_vec);
  R_xlen_t new_len = XLENGTH(new_vals);
  SEXP out;

  if (old_type == STRSXP && new_type == STRSXP) {
    out = PROTECT(Rf_allocVector(STRSXP, old_len + new_len));
    for (R_xlen_t pos = 0; pos < old_len; pos++) {
      SET_STRING_ELT(out, pos, STRING_ELT(old_vec, pos));
    }
    for (R_xlen_t pos = 0; pos < new_len; pos++) {
      SET_STRING_ELT(out, old_len + pos, STRING_ELT(new_vals, pos));
    }
  } else if (old_type == REALSXP && new_type == REALSXP) {
    out = PROTECT(Rf_allocVector(REALSXP, old_len + new_len));
    memcpy(REAL(out), REAL(old_vec), old_len * sizeof(double));
    memcpy(REAL(out) + old_len, REAL(new_vals), new_len * sizeof(double));
  } else if (old_type == REALSXP && (new_type == INTSXP || new_type == LGLSXP)) {
    out = PROTECT(Rf_allocVector(REALSXP, old_len + new_len));
    double *out_ptr = REAL(out);
    const int *new_ptr = (new_type == INTSXP) ? INTEGER(new_vals) : LOGICAL(new_vals);
    memcpy(out_ptr, REAL(old_vec), old_len * sizeof(double));
    for (R_xlen_t pos = 0; pos < new_len; pos++) {
      out_ptr[old_len + pos] =
        (new_ptr[pos] == NA_INTEGER) ? NA_REAL : (double) new_ptr[pos];
    }
  } else if ((old_type == INTSXP || old_type == LGLSXP) &&
             (new_type == old_type || new_type == LGLSXP)) {
    /* integer and logical vectors both store int */
    out = PROTECT(Rf_allocVector(old_type, old_len + new_len));
    int *out_ptr = (old_type == INTSXP) ? INTEGER(out) : LOGICAL(out);
    const int *old_ptr = (old_type == INTSXP) ? INTEGER(old_vec) : LOGICAL(old_vec);
    const int *new_ptr = (new_type == INTSXP) ? INTEGER(new_vals) : LOGICAL(new_vals);
    memcpy(out_ptr, old_ptr, old_len * sizeof(int));
    memcpy(out_ptr + old_len, new_ptr, new_len * sizeof(int));
  } else {
    Rf_error("c() does not keep the type of old_vec for this pair of types");
  }
  UNPROTECT(1);
  return out;
}
