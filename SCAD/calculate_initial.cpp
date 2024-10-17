#include <Rcpp.h>
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
Eigen::VectorXd calculate_initial(Eigen::MatrixXd X, Eigen::VectorXd Y) {
  int n = X.rows();
  int p = X.cols();
  
  Eigen::VectorXd ans = Eigen::VectorXd::Zero(p);
  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if(i != j) {
        ans += L_prime(Y[i] - Y[j]) * (X.row(i) - X.row(j));
      }
    }
  }
  
  ans /= (n * (n - 1));
  
  return ans;
}
