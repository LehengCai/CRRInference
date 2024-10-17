#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd calculate_ck(Eigen::MatrixXd X) {
  int n = X.rows();
  int p = X.cols();
  
  Eigen::VectorXd c_k = Eigen::VectorXd::Zero(p);
  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      if(i != j) {
        c_k += (X.row(i) - X.row(j)).array().square().matrix();
      }
    }
  }
  
  c_k /= (n * (n - 1));
  
  return c_k;
}
