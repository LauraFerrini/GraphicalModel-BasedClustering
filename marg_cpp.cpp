#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double marg_S_test(IntegerVector S,
              DataFrame N,
              double a,
              IntegerVector I_cal) {
  
  // If S is empty, return 0
  if (S.size() == 0) {
    return 0.0;
  }
  
  // Get number of levels for the current set of variables S
  int prod_levels = 1;
  IntegerVector S0 = clone(S) - 1; // Convert to 0-based for C++
  for (int i = 0; i < S.size(); ++i) {
    prod_levels *= I_cal[S[i] - 1];
  }
  
  // Create a map from combinations in S to their total frequencies
  std::map<std::vector<int>, double> count_map;
  
  IntegerVector freq = N["freq"];
  int n_rows = freq.size();
  
  // Extract the subset columns from N
  std::vector< IntegerVector > subset_cols;
  for (int i = 0; i < S.size(); ++i) {
    subset_cols.push_back(N[S[i] - 1]);  // Column S[i] (1-based)
  }
  
  for (int i = 0; i < n_rows; ++i) {
    std::vector<int> key(S.size());
    for (int j = 0; j < S.size(); ++j) {
      key[j] = subset_cols[j][i];
    }
    count_map[key] += freq[i];
  }
  
  // Build full count vector (observed + zeroes)
  std::vector<double> N_S;
  for (auto& pair : count_map) {
    N_S.push_back(pair.second);
  }
  
  int n_obs = N_S.size();
  int n_zero = prod_levels - n_obs;
  for (int i = 0; i < n_zero; ++i) {
    N_S.push_back(0.0);
  }
  
  double prior = a / static_cast<double>(prod_levels);
  double total = 0.0;
  for (size_t i = 0; i < N_S.size(); ++i) {
    total += R::lgammafn(prior + N_S[i]) - R::lgammafn(prior);
  }
  
  double result = R::lgammafn(a) - R::lgammafn(a + std::accumulate(N_S.begin(), N_S.end(), 0.0)) + total;
  
  return result;
}



// [[Rcpp::export]]
double logprior_ratio(IntegerVector G_star,
                      IntegerVector G,
                      double a_pi,
                      double b_pi,
                      int q) {
  
  int sum_G_star = std::accumulate(G_star.begin(), G_star.end(), 0);
  int sum_G = std::accumulate(G.begin(), G.end(), 0);
  int total_edges = q * (q - 1) / 2;
  
  double logprior_star = R::lgammafn(sum_G_star / 2.0 + a_pi) +
    R::lgammafn(total_edges - sum_G_star / 2.0 + b_pi);
  
  double logprior = R::lgammafn(sum_G / 2.0 + a_pi) +
    R::lgammafn(total_edges - sum_G / 2.0 + b_pi);
  
  return logprior_star - logprior;
}
