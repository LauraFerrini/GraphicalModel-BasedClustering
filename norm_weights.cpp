#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector normalize_weights_test(NumericVector prob) {
  int n = prob.size();
  NumericVector weight(n);
  
  // Step 1: collect values that are not -Inf
  std::vector<double> finite_vals;
  for (int i = 0; i < n; ++i) {
    if (prob[i] != R_NegInf) {
      finite_vals.push_back(prob[i]);
    }
  }
  
  // Step 2: compute mean of finite values (like in R code)
  double const_term = 0.0;
  if (!finite_vals.empty()) {
    const_term = std::accumulate(finite_vals.begin(), finite_vals.end(), 0.0) / finite_vals.size();
  }
  
  // Step 3: compute weights
  double denom = 0.0;
  for (int i = 0; i < n; ++i) {
    if (prob[i] != R_NegInf) {
      double temp = std::exp(prob[i] + const_term);
      weight[i] = temp;
      denom += temp;
    }
  }
  
  // Step 4: normalize
  if (denom > 0) {
    for (int i = 0; i < n; ++i) {
      weight[i] = (prob[i] == R_NegInf) ? 0.0 : weight[i] / denom;
    }
  }
  
  return weight;
}