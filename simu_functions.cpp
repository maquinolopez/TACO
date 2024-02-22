#include <Rcpp.h>
using namespace Rcpp;


// Simple linear interpolation function
double linearInterpolate(double x, double x0, double x1, double y0, double y1) {
  if (x1 == x0) {
    return (y0 + y1) / 2; // Avoid division by zero; return the average or some other sensible value
  }
  return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

// age_depth function using linear interpolation
// Assumes bot_depth and bot_ages are accessible and have size of at least 2
double age_depth(double x, NumericVector bot_depth, NumericVector bot_ages) {
  // Example implementation; actual implementation may need to handle more cases
  return linearInterpolate(x, bot_depth[0], bot_depth[1], bot_ages[0], bot_ages[1]);
}


// [[Rcpp::export]]
NumericVector w(NumericVector param, NumericVector depths) { // Pondering function (it is a logistic regression with a fix)
  double k = param[3] / param[4]; // Adjusted index for C++ (0-based indexing)
  double v = param[4];            // Adjusted index for C++ (0-based indexing)
  NumericVector ws(depths.size());
  
  for(int i = 0; i < depths.size(); i++) {
    ws[i] = 1 / (1 + exp(-k * (depths[i] - v)));
  }
  
  return ws;
}



// [[Rcpp::export]]
double M_t(double m0, double alpha, double t_x, double t_x_delta) {
  double Mi = (-m0 / alpha) * (exp(-alpha * t_x) - exp(-alpha * t_x_delta));
  return Mi;
}


// [[Rcpp::export]]
DataFrame m0s(NumericVector param, NumericVector bot_depth, NumericVector bot_ages) {
  int len_data = bot_depth.size(); // This replaces the len_data in R, ensuring it matches the depth array
  NumericVector T_a(len_data), T_c(len_data), m0_a(len_data), m0_c(len_data);
  
  for(int i = 0; i < len_data; i++) {
    T_a[i] = age_depth(param[4] - param[3], bot_depth, bot_ages); // Top age
    T_c[i] = age_depth(param[4] + param[3], bot_depth, bot_ages); // Bottom age
    
    // Assuming M_t can be vectorized or adjusted accordingly
    m0_a[i] = M_t(param[0], param[1], T_a[i], 0);
    m0_c[i] = M_t(m0_a[i], param[2], T_c[i] - T_a[i], 0);
  }
  
  // Creating a DataFrame to return
  return DataFrame::create(Named("m0") = rep(param[0], len_data), 
                           Named("m0a") = m0_a, 
                           Named("m0c") = m0_c);
}


// [[Rcpp::export]]
NumericVector simu(NumericVector param, NumericVector bot_depth, NumericVector bot_ages, NumericVector top_ages) {
  double T_a = param[4] - param[3]; // Adjusted for 0-based indexing
  double T_c = param[4] + param[3]; // Adjusted for 0-based indexing
  
  NumericVector ws = w(param); // Assuming w is adapted for C++ and vectorized appropriately
  
  // Determine indices for each phase
  IntegerVector indx_phase1, indx_phase2, indx_phase3;
  for(int i = 0; i < bot_depth.size(); i++) {
    if (bot_depth[i] <= T_a) indx_phase1.push_back(i);
    else if (bot_depth[i] > T_a && bot_depth[i] <= T_c) indx_phase2.push_back(i);
    else if (bot_depth[i] > T_c) indx_phase3.push_back(i);
  }
  
  // Calculate M_0s
  DataFrame M_0s = m0s(param); // Assuming m0s adapted for C++ and returns a DataFrame
  
  NumericVector mi_1, mi_2, mi_3; // To store results for each phase
  
  // Phase 1 calculations (simplified for clarity)
  for(int i : indx_phase1) {
    mi_1.push_back(M_t(M_0s["m0"][i], param[1], bot_ages[i], top_ages[i]));
  }
  
  // Phase 2 and 3 calculations would follow a similar pattern, adjusted for logic and parameters
  
  // Combine and return results
  NumericVector z = concat(mi_1, mi_2, mi_3); // You'll need to implement concat or use an alternative approach
  return z;
}

