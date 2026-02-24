#include <Rcpp.h>
using namespace Rcpp;

// remove_node_attr_cpp
//
// Efficiently removes elements at specified positions from every vector in
// a named list of node attributes. Replaces the R pattern:
//   lapply(attr_list, function(x) x[-posit_ids])
//
// By pre-computing a boolean keep-mask once and applying it across all
// attribute vectors, we avoid repeated interpretation of the negative-index
// vector for each attribute.
//
// Handles integer, double, character, and logical vectors (the types used
// by EpiModel node attributes).
//
// @param attr_list Named list of attribute vectors (all same length)
// @param posit_ids Integer vector of 1-indexed positions to remove
// @return A new named list with the specified positions removed from each vector
//
// [[Rcpp::export]]
List remove_node_attr_cpp(List attr_list, IntegerVector posit_ids) {

  int n_attr = attr_list.size();

  if (posit_ids.size() == 0 || n_attr == 0) {
    return clone(attr_list);
  }

  // Get the length from the first attribute vector
  SEXP first = attr_list[0];
  int n = Rf_length(first);

  // Build boolean keep mask (pre-computed once, applied to all attributes)
  std::vector<bool> keep(n, true);
  for (int i = 0; i < posit_ids.size(); i++) {
    int idx = posit_ids[i] - 1; // convert from R's 1-indexing
    if (idx >= 0 && idx < n) {
      keep[idx] = false;
    }
  }

  int n_keep = 0;
  for (int i = 0; i < n; i++) {
    if (keep[i]) n_keep++;
  }

  List result(n_attr);
  result.names() = attr_list.names();

  for (int a = 0; a < n_attr; a++) {
    SEXP elem = attr_list[a];
    switch (TYPEOF(elem)) {
      case INTSXP: {
        IntegerVector v(elem);
        IntegerVector out(n_keep);
        int j = 0;
        for (int i = 0; i < n; i++) {
          if (keep[i]) out[j++] = v[i];
        }
        result[a] = out;
        break;
      }
      case REALSXP: {
        NumericVector v(elem);
        NumericVector out(n_keep);
        int j = 0;
        for (int i = 0; i < n; i++) {
          if (keep[i]) out[j++] = v[i];
        }
        result[a] = out;
        break;
      }
      case STRSXP: {
        CharacterVector v(elem);
        CharacterVector out(n_keep);
        int j = 0;
        for (int i = 0; i < n; i++) {
          if (keep[i]) out[j++] = v[i];
        }
        result[a] = out;
        break;
      }
      case LGLSXP: {
        LogicalVector v(elem);
        LogicalVector out(n_keep);
        int j = 0;
        for (int i = 0; i < n; i++) {
          if (keep[i]) out[j++] = v[i];
        }
        result[a] = out;
        break;
      }
      default: {
        // Fallback for any other type: build a LogicalVector keep-mask
        // and use R's subsetting operator
        LogicalVector keep_lgl(n);
        for (int i = 0; i < n; i++) {
          keep_lgl[i] = keep[i];
        }
        Environment base("package:base");
        Function subset = base["["];
        result[a] = subset(elem, keep_lgl);
        break;
      }
    }
  }

  return result;
}
