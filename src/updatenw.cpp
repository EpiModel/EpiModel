#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// assumes x and y are sorted in increasing order

// [[Rcpp::export]]
IntegerVector shiftVec(IntegerVector x, IntegerVector y) {

  int x_len = x.size();
  int y_len = y.size();

  IntegerVector z(x_len);

  int i = 0;

  for(int j = 0; j < y_len; j++) {
    while(i < x_len && x[i] < y[j]) {
      z[i] = x[i] - j;
      i++;
    }
  }

  // handle any trailing part of x after getting to the end of y
  while(i < x_len) {
    z[i] = x[i] - y_len;
    i++;
  }

  return(z);
}

// delete_vertices_cpp
//
// Consolidated C++ implementation that replaces the R pipeline of:
//   delete_edges() + 4x order() + shiftVec() + 4x order() unsort
//
// Given an edgelist matrix (2-col integer matrix) and a sorted vector of
// vertex IDs to delete:
//   1. Removes any edges involving the deleted vertex IDs
//   2. Shifts remaining vertex IDs downward to fill gaps
//
// Uses a precomputed shift-amount lookup table instead of sorting, reducing
// complexity from O(E log E) to O(N + E).
//
// @param el Two-column integer matrix of edges (1-indexed node IDs)
// @param vid Sorted integer vector of vertex IDs to delete (1-indexed)
// @param n Total network size before deletion
// @return A two-column integer matrix with edges updated
//
// [[Rcpp::export]]
IntegerMatrix delete_vertices_cpp(IntegerMatrix el, IntegerVector vid, int n) {

  int nEdges = el.nrow();
  int nDel = vid.size();

  if (nDel == 0) {
    // Nothing to delete -- return a copy
    IntegerMatrix out = clone(el);
    return out;
  }

  if (nEdges == 0) {
    // No edges -- return empty 0x2 matrix
    IntegerMatrix out(0, 2);
    return out;
  }

  // Build a boolean set for O(1) vertex deletion lookup
  // and a shift table: shift[v] = number of deleted IDs <= v
  // Both are 1-indexed (index 0 unused)
  std::vector<bool> is_deleted(n + 1, false);
  for (int i = 0; i < nDel; i++) {
    is_deleted[vid[i]] = true;
  }

  // Precompute shift amounts: shift[v] = how many deleted IDs are < v
  // This replaces the 4x sort + shiftVec approach
  std::vector<int> shift(n + 1, 0);
  int cumulative = 0;
  for (int v = 1; v <= n; v++) {
    if (is_deleted[v]) {
      cumulative++;
    }
    shift[v] = cumulative;
  }

  // First pass: count surviving edges
  int nSurvive = 0;
  for (int i = 0; i < nEdges; i++) {
    int n1 = el(i, 0);
    int n2 = el(i, 1);
    if (!is_deleted[n1] && !is_deleted[n2]) {
      nSurvive++;
    }
  }

  // Second pass: compact edges and apply shift
  IntegerMatrix out(nSurvive, 2);
  int idx = 0;
  for (int i = 0; i < nEdges; i++) {
    int n1 = el(i, 0);
    int n2 = el(i, 1);
    if (!is_deleted[n1] && !is_deleted[n2]) {
      out(idx, 0) = n1 - shift[n1];
      out(idx, 1) = n2 - shift[n2];
      idx++;
    }
  }

  return out;
}
