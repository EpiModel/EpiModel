#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP shiftVec(SEXP x, SEXP y);
SEXP compact_attr_list(SEXP x, SEXP drop);
SEXP append_attr_c(SEXP old, SEXP nw);

static const R_CallMethodDef CallEntries[] = {
  {"shiftVec", (DL_FUNC) &shiftVec, 2},
  {"compact_attr_list", (DL_FUNC) &compact_attr_list, 2},
  {"append_attr_c", (DL_FUNC) &append_attr_c, 2},
  {NULL, NULL, 0}
};

void R_init_EpiModel(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  /* ergm looks up the change statistics in changestats.users.c by symbol
     name, so dynamic symbol lookup must stay on */
  R_useDynamicSymbols(dll, TRUE);
}
