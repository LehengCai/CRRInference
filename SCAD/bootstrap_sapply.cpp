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
NumericMatrix bootstrap_sapply(int B, Eigen::VectorXd epsilon_hat, Eigen::MatrixXd X, Eigen::MatrixXd W_hat, Eigen::MatrixXd S_hat_sqrt_inverse) {
  int n = X.rows();
  int p = X.cols();
  
  NumericMatrix bootstraps(3, B);
  
  for(int b = 0; b < B; ++b) {
	// Generate standard normal random numbers using rnorm
    NumericVector z = rnorm(n);
    // Convert NumericVector to Eigen VectorXd
    Map<VectorXd> Z(as<Map<VectorXd> >(z));//Eigen::VectorXd Z = Eigen::VectorXd::Random(n);
    Eigen::VectorXd ans = Eigen::VectorXd::Zero(p);
    
    for(int i = 0; i < n; ++i) {
      for(int j = 0; j < n; ++j) {
        if(i != j) {
          ans += L_prime(epsilon_hat[i] - epsilon_hat[j]) * (X.row(i) - X.row(j)) * (Z[i] + Z[j]);
        }
      }
    }
    
    ans = S_hat_sqrt_inverse * W_hat * ans / (n * (n - 1));
    
    double max1 = ans.head(5).cwiseAbs().maxCoeff();
    double max2 = ans.head(floor(p / 5)).cwiseAbs().maxCoeff();
    double max3 = ans.cwiseAbs().maxCoeff();
    
    bootstraps(0, b) = max1;
    bootstraps(1, b) = max2;
    bootstraps(2, b) = max3;
  }
  
  return bootstraps;
}
