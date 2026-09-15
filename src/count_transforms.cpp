#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cpp_rel_abundance(NumericMatrix counts) {
  int nrow = counts.nrow();
  int ncol = counts.ncol();
  
  NumericMatrix result(nrow, ncol);
  
  for (int j = 0; j < ncol; j++) {
    double col_sum = 0;
    
    // Calculate column sum
    for (int i = 0; i < nrow; i++) {
      col_sum += counts(i, j);
    }
    
    // Calculate relative abundance
    if (col_sum > 0) {
      for (int i = 0; i < nrow; i++) {
        result(i, j) = counts(i, j) / col_sum;
      }
    } else {
      // Handle zero sums (set all to 0)
      for (int i = 0; i < nrow; i++) {
        result(i, j) = 0.0;
      }
    }
  }
  
  // Set row and column names
  result.attr("dimnames") = counts.attr("dimnames");
  
  return result;
}

// [[Rcpp::export]]
NumericMatrix cpp_clr_transform(NumericMatrix counts, double pseudocount = 1.0) {
  int nrow = counts.nrow();
  int ncol = counts.ncol();
  
  NumericMatrix result(nrow, ncol);
  
  for (int j = 0; j < ncol; j++) {
    // Apply pseudocount and take log
    NumericVector log_counts(nrow);
    double sum_log = 0.0;
    
    for (int i = 0; i < nrow; i++) {
      log_counts[i] = log(counts(i, j) + pseudocount);
      sum_log += log_counts[i];
    }
    
    // Calculate geometric mean
    double geom_mean = sum_log / nrow;
    
    // Calculate CLR transform
    for (int i = 0; i < nrow; i++) {
      result(i, j) = log_counts[i] - geom_mean;
    }
  }
  
  // Set row and column names
  result.attr("dimnames") = counts.attr("dimnames");
  
  return result;
}

// [[Rcpp::export]]
NumericMatrix cpp_log_transform(NumericMatrix counts, double pseudocount = 1.0) {
  int nrow = counts.nrow();
  int ncol = counts.ncol();
  
  NumericMatrix result(nrow, ncol);
  
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      result(i, j) = log(counts(i, j) + pseudocount);
    }
  }
  
  // Set row and column names
  result.attr("dimnames") = counts.attr("dimnames");
  
  return result;
}

// [[Rcpp::export]]
LogicalMatrix cpp_presence_absence(NumericMatrix counts) {
  int nrow = counts.nrow();
  int ncol = counts.ncol();
  
  LogicalMatrix result(nrow, ncol);
  
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      result(i, j) = counts(i, j) > 0;
    }
  }
  
  // Set row and column names
  result.attr("dimnames") = counts.attr("dimnames");
  
  return result;
}