#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// update_cuml_edgelist_cpp
//
// Purpose-built hash join that replaces the dplyr::full_join() used in
// update_cumulative_edgelist(). Given a cumulative active edgelist and
// the current edgelist (both identified by unique ID pairs), this function
// categorizes edges into three groups in a single pass:
//
//   1. Continuing edges: present in both cumulative and current
//   2. New edges: present in current but not cumulative (assigned start = at)
//   3. Terminated edges: present in cumulative but not current (assigned stop = at-1)
//
// Uses an unordered_map with a pair hash for O(1) lookup per edge,
// reducing complexity from dplyr's general-purpose join to a specialized
// O(E_cuml + E_cur) operation.
//
// @param cuml_head,cuml_tail,cuml_start Vectors from the cumulative active
//   edgelist (unique IDs stored as numeric/double)
// @param cur_head,cur_tail Vectors from the current edgelist (unique IDs)
// @param at Current time step (integer)
// @return A list with 7 elements: active_head, active_tail, active_start
//   (continuing + new edges), and term_head, term_tail, term_start, term_stop
//   (terminated edges)
//

// Hash function for (int64, int64) pairs used as edge keys
struct PairHash {
  std::size_t operator()(const std::pair<int64_t, int64_t>& p) const {
    // Combine hashes using a multiplicative mixing approach
    std::size_t h1 = std::hash<int64_t>{}(p.first);
    std::size_t h2 = std::hash<int64_t>{}(p.second);
    return h1 ^ (h2 * 2654435761ULL);
  }
};

// [[Rcpp::export]]
List update_cuml_edgelist_cpp(NumericVector cuml_head, NumericVector cuml_tail,
                              NumericVector cuml_start,
                              NumericVector cur_head, NumericVector cur_tail,
                              int at) {

  int n_cuml = cuml_head.size();
  int n_cur = cur_head.size();

  // Build hash map from cumulative edges: (head, tail) -> index
  // Use int64_t keys since unique IDs are integer-valued doubles
  typedef std::pair<int64_t, int64_t> EdgeKey;
  std::unordered_map<EdgeKey, int, PairHash> cuml_map;
  cuml_map.reserve(n_cuml);

  for (int i = 0; i < n_cuml; i++) {
    EdgeKey key(static_cast<int64_t>(cuml_head[i]),
                static_cast<int64_t>(cuml_tail[i]));
    cuml_map[key] = i;
  }

  // Track which cumulative edges are still active (matched by current edges)
  std::vector<bool> cuml_matched(n_cuml, false);

  // Track which current edges are new (not in cumulative)
  std::vector<bool> cur_is_new(n_cur, false);

  for (int j = 0; j < n_cur; j++) {
    EdgeKey key(static_cast<int64_t>(cur_head[j]),
                static_cast<int64_t>(cur_tail[j]));
    auto it = cuml_map.find(key);
    if (it != cuml_map.end()) {
      cuml_matched[it->second] = true;
    } else {
      cur_is_new[j] = true;
    }
  }

  // Count results for pre-allocation
  int n_continuing = 0;
  int n_new = 0;
  int n_terminated = 0;

  for (int i = 0; i < n_cuml; i++) {
    if (cuml_matched[i]) n_continuing++;
    else n_terminated++;
  }
  for (int j = 0; j < n_cur; j++) {
    if (cur_is_new[j]) n_new++;
  }

  int n_active = n_continuing + n_new;

  // Build active output: continuing edges first, then new edges
  NumericVector act_head(n_active);
  NumericVector act_tail(n_active);
  NumericVector act_start(n_active);

  int idx = 0;

  // Continuing edges (preserve original start time)
  for (int i = 0; i < n_cuml; i++) {
    if (cuml_matched[i]) {
      act_head[idx] = cuml_head[i];
      act_tail[idx] = cuml_tail[i];
      act_start[idx] = cuml_start[i];
      idx++;
    }
  }

  // New edges (start = current time step)
  for (int j = 0; j < n_cur; j++) {
    if (cur_is_new[j]) {
      act_head[idx] = cur_head[j];
      act_tail[idx] = cur_tail[j];
      act_start[idx] = static_cast<double>(at);
      idx++;
    }
  }

  // Build terminated output (stop = at - 1)
  NumericVector term_head(n_terminated);
  NumericVector term_tail(n_terminated);
  NumericVector term_start(n_terminated);
  NumericVector term_stop(n_terminated);

  idx = 0;
  for (int i = 0; i < n_cuml; i++) {
    if (!cuml_matched[i]) {
      term_head[idx] = cuml_head[i];
      term_tail[idx] = cuml_tail[i];
      term_start[idx] = cuml_start[i];
      term_stop[idx] = static_cast<double>(at - 1);
      idx++;
    }
  }

  return List::create(
    Named("active_head") = act_head,
    Named("active_tail") = act_tail,
    Named("active_start") = act_start,
    Named("term_head") = term_head,
    Named("term_tail") = term_tail,
    Named("term_start") = term_start,
    Named("term_stop") = term_stop
  );
}
