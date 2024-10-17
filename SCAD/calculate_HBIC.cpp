#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// Define the function L
double L(double u) {
  if(u >= -1 && u <= 1) {
    return 3.0/4.0 * u * u - u * u * u * u / 8.0 + 3.0/8.0;
  } else if (u > 1) {
    return u;
  } else {
    return -u;
  }
}

// Define the function L_h
double L_h(double u, double h = 1) {
  return h * L(u / h);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericVector calculate_HBIC(Eigen::MatrixXd solutions, Eigen::MatrixXd X, Eigen::VectorXd Y) {
  int n = X.rows();
  int p = X.cols();
  int num_solutions = solutions.rows();
  
  NumericVector HBIC(num_solutions);
  
  for(int k = 0; k < num_solutions; ++k) {
    double ans = 0;
    Eigen::VectorXd b = solutions.row(k);
    
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j) {
        if(i != j) {
          ans += L_h(Y[i] - Y[j] - (X.row(i) - X.row(j)).dot(b));
        }
      }
    }
    
    ans /= (n * (n - 1));
    HBIC[k] = log(ans) + log(log(n)) * log(p) / n * (b.array() != 0).count();
  }
  
  return HBIC;
}
