#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// discord_edgelist_cpp
//
// Single-pass C++ implementation of discord_edgelist(). Given an edgelist
// matrix, node status vector, active vector, and a set of infectious statuses,
// identifies all discordant (susceptible-infected) pairs. The edgelist rows are
// processed in random order (Fisher-Yates shuffle). Both endpoints must be
// active (== 1).
//
// Returns a List with three integer vectors: at, sus, inf. These can be
// assembled into a data.frame by the R wrapper.
//
// [[Rcpp::export]]
List discord_edgelist_cpp(IntegerMatrix el,
                          CharacterVector status,
                          IntegerVector active,
                          CharacterVector infstat,
                          int at_val) {

  int nEdges = el.nrow();

  if (nEdges == 0) {
    return List::create(
      Named("at")  = IntegerVector(0),
      Named("sus") = IntegerVector(0),
      Named("inf") = IntegerVector(0)
    );
  }

  // Build a set of infectious status values for O(1) lookup
  std::unordered_set<std::string> infset;
  for (int i = 0; i < infstat.size(); i++) {
    infset.insert(std::string(infstat[i]));
  }

  // Create a shuffled index array (Fisher-Yates)
  std::vector<int> idx(nEdges);
  for (int i = 0; i < nEdges; i++) idx[i] = i;
  for (int i = nEdges - 1; i > 0; i--) {
    int j = (int)(R::runif(0.0, 1.0) * (i + 1));
    if (j > i) j = i;  // guard against edge case
    std::swap(idx[i], idx[j]);
  }

  // Pre-allocate output vectors (worst case: every edge is discordant)
  std::vector<int> sus_out;
  std::vector<int> inf_out;
  sus_out.reserve(nEdges);
  inf_out.reserve(nEdges);

  // Single pass over shuffled edges
  for (int k = 0; k < nEdges; k++) {
    int row = idx[k];
    int n1 = el(row, 0);  // 1-indexed node IDs
    int n2 = el(row, 1);

    // Check both nodes are active (R vectors are 0-indexed in C++)
    if (active[n1 - 1] != 1 || active[n2 - 1] != 1) continue;

    std::string s1(status[n1 - 1]);
    std::string s2(status[n2 - 1]);

    bool s1_is_sus = (s1 == "s");
    bool s2_is_sus = (s2 == "s");
    bool s1_is_inf = (infset.count(s1) > 0);
    bool s2_is_inf = (infset.count(s2) > 0);

    // S-I pair: node1 susceptible, node2 infected
    if (s1_is_sus && s2_is_inf) {
      sus_out.push_back(n1);
      inf_out.push_back(n2);
    }
    // I-S pair: node1 infected, node2 susceptible (swap to sus, inf order)
    else if (s2_is_sus && s1_is_inf) {
      sus_out.push_back(n2);
      inf_out.push_back(n1);
    }
  }

  int nPairs = sus_out.size();
  if (nPairs == 0) {
    return List::create(
      Named("at")  = IntegerVector(0),
      Named("sus") = IntegerVector(0),
      Named("inf") = IntegerVector(0)
    );
  }

  IntegerVector at_out(nPairs, at_val);
  IntegerVector sus_vec(sus_out.begin(), sus_out.end());
  IntegerVector inf_vec(inf_out.begin(), inf_out.end());

  return List::create(
    Named("at")  = at_out,
    Named("sus") = sus_vec,
    Named("inf") = inf_vec
  );
}
