#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// Define the function L_prime
double L_prime(double u) {
  if(u >= -1 && u <= 1) {
    return 3.0/2.0*u - pow(u, 3)/2;
  } else if (u > 1) {
    return 1.0;
  } else {
    return -1.0;
  }
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd calculate_v_and_a(Eigen::MatrixXd X, Eigen::VectorXd Y, Eigen::VectorXd beta_tilde) {
  int n = X.rows();
  int p = X.cols();
  
  Eigen::MatrixXd v(n, n);
  Eigen::VectorXd a(p);
  
  // Calculate v matrix
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      v(i, j) = Y[i] - Y[j] - (X.row(i) - X.row(j)).dot(beta_tilde);
    }
  }
  
  // Initialize a to zero
  a.setZero();
  
  // Calculate a vector
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if(i != j) {
        a += L_prime(v(i, j)) * (X.row(i) - X.row(j)).transpose();
      }
    }
  }
  
  // Normalize a
  a = -a / (n * (n - 1));
  
  return a;
}
