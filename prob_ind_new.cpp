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
IntegerVector compareToRow(NumericMatrix mat, NumericVector x_subset) {
  int n = mat.nrow();
  int p = mat.ncol();
  IntegerVector match(n, 1);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      if (mat(i, j) != x_subset[j]) {
        match[i] = 0;
        break;
      }
    }
  }
  
  return match;
}

// [[Rcpp::export]]
double prob_ik_nonempty_test(NumericMatrix Xk,
                        List cl_list,
                        List se_list,
                        bool member,
                        NumericVector x,
                        double a,
                        IntegerVector I_cal) {
  
  List se_list_filtered;
  for (int i = 0; i < se_list.size(); ++i) {
    if (Rf_length(se_list[i]) > 0) {
      se_list_filtered.push_back(se_list[i]);
    }
  }
  
  int n = Xk.nrow();
  double log_n_member = std::log(n - (int)member);
  
  double cl_term = 0.0;
  for (int i = 0; i < cl_list.size(); ++i) {
    IntegerVector j = cl_list[i];
    NumericMatrix sub_Xk(n, j.size());
    NumericVector xj(j.size());
    
    for (int col = 0; col < j.size(); ++col) {
      int idx = j[col] - 1;
      xj[col] = x[idx];
      for (int row = 0; row < n; ++row) {
        sub_Xk(row, col) = Xk(row, idx);
      }
    }
    
    IntegerVector match = compareToRow(sub_Xk, xj);
    int count = std::accumulate(match.begin(), match.end(), 0);
    
    double denom = 1.0;
    for (int d = 0; d < j.size(); ++d) {
      denom *= I_cal[j[d] - 1];
    }
    
    cl_term += std::log(a / denom + count - (int)member);
  }
  
  double se_term = 0.0;
  for (int i = 0; i < se_list_filtered.size(); ++i) {
    IntegerVector j = se_list_filtered[i];
    NumericMatrix sub_Xk(n, j.size());
    NumericVector xj(j.size());
    
    for (int col = 0; col < j.size(); ++col) {
      int idx = j[col] - 1;
      xj[col] = x[idx];
      for (int row = 0; row < n; ++row) {
        sub_Xk(row, col) = Xk(row, idx);
      }
    }
    
    IntegerVector match = compareToRow(sub_Xk, xj);
    int count = std::accumulate(match.begin(), match.end(), 0);
    
    double denom = 1.0;
    for (int d = 0; d < j.size(); ++d) {
      denom *= I_cal[j[d] - 1];
    }
    
    se_term += std::log(a / denom + count - (int)member);
  }
  
  return log_n_member + (se_list_filtered.size() - cl_list.size()) * std::log(a + n - (int)member) + cl_term - se_term;
}




// [[Rcpp::export]]
double prob_ik_empty_test(List cl_list,
                     List se_list,
                     double a,
                     double alpha0,
                     IntegerVector I_cal) {
  
  // Filter out empty elements from se_list
  List se_list_filtered;
  for (int i = 0; i < se_list.size(); ++i) {
    if (Rf_length(se_list[i]) > 0) {
      se_list_filtered.push_back(se_list[i]);
    }
  }
  
  double cl_term = 0.0;
  for (int i = 0; i < cl_list.size(); ++i) {
    IntegerVector j = cl_list[i];
    double prod_I = 1.0;
    for (int d = 0; d < j.size(); ++d) {
      prod_I *= I_cal[j[d] - 1];  // Convert to 0-based index
    }
    cl_term += std::log(1.0 / prod_I);
  }
  
  double se_term = 0.0;
  for (int i = 0; i < se_list_filtered.size(); ++i) {
    IntegerVector j = se_list_filtered[i];
    double prod_I = 1.0;
    for (int d = 0; d < j.size(); ++d) {
      prod_I *= I_cal[j[d] - 1];
    }
    se_term += std::log(1.0 / prod_I);
  }
  
  double result = std::log(alpha0) + (se_list_filtered.size() - cl_list.size()) * std::log(a) + cl_term - se_term;
  return result;
}



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



