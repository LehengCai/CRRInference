library(flare)
library(Rcpp)
library(RcppEigen)
sourceCpp('bootstrap_sapply.cpp')

#inference for the first five coefficients, the first [p/5] coefficients, and all the coefficients.

inference =  function(X,Y,beta_hat,beta_true,alpha=0.05,B=500){
  n = nrow(X)
  p = ncol(X)
  Sigma_hat_inverse = sugm(X, verbose = F,nlambda = 2, lambda.min.ratio = 0.25)#0.2
  Sigma_hat_inverse = Sigma_hat_inverse$icov[[Sigma_hat_inverse$nlambda]]
  
  epsilon_hat = Y-X%*%beta_hat
  
  denominator = 0
  ans = rep(0,p)
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        ans = ans + L_prime(epsilon_hat[i] - epsilon_hat[j])*(X[i,]-X[j,])
        denominator = denominator+K_h(epsilon_hat[i]-epsilon_hat[j])
      }
    }
  }
  ans = ans/(n*(n-1))
  denominator = 4*denominator / (n*(n-1))
  W_hat = Sigma_hat_inverse/denominator
  S_hat_sqrt_diag  =  sqrt(diag(W_hat))
  S_hat_sqrt_inverse = matrix(0,p,p)
  diag(S_hat_sqrt_inverse) =   1/sqrt(diag(W_hat))  
  
  debias = beta_hat + W_hat%*% ans
  
  bootstraps = bootstrap_sapply(B=B,epsilon_hat,X,W_hat,S_hat_sqrt_inverse)
  Q = apply(bootstraps,1,quantile,1-alpha)
  
  # bootstraps = sapply(1:B,function(xxx,epsilon_hat,X,W_hat,S_hat_sqrt_inverse){
  #   n = nrow(X)
  #   p = ncol(X)
  #   Z = rnorm(n)
  #   ans = rep(0,p)
  #   for(i in 1:n){
  #     for(j in 1:n){
  #       if(i != j){
  #         ans = ans + L_prime(epsilon_hat[i] - epsilon_hat[j])*
  #           (X[i,]-X[j,])*(Z[i]+Z[j])
  #       }
  #     }
  #   }
  #   ans = S_hat_sqrt_inverse%*%W_hat%*%ans/(n*(n-1))
  #   return(c(max(abs(ans[1:5])), max(abs(ans[1:floor(p/5)])), max(abs(ans))))
  # },epsilon_hat,X,W_hat,S_hat_sqrt_inverse)
  # Q = apply(bootstraps,1,quantile,1-alpha)

  #beta_true = c(rep(sqrt(3),3),rep(0,p-3))
  return(c(all((debias[1:5]-S_hat_sqrt_diag[1:5]*Q[1] < beta_true[1:5]) &
                 (debias[1:5]+S_hat_sqrt_diag[1:5]*Q[1] > beta_true[1:5])),
           mean(S_hat_sqrt_diag[1:5]*Q[1]*2),
           all((debias[1:floor(p/5)]-S_hat_sqrt_diag[1:floor(p/5)]*Q[2] < beta_true[1:floor(p/5)]) &
                 (debias[1:floor(p/5)]+S_hat_sqrt_diag[1:floor(p/5)]*Q[2] > beta_true[1:floor(p/5)])),
           mean(S_hat_sqrt_diag[1:floor(p/5)]*Q[2]*2),
           all((debias - S_hat_sqrt_diag*Q[3] < beta_true ) &
                 (debias + S_hat_sqrt_diag*Q[3] > beta_true )),
           mean(S_hat_sqrt_diag*Q[3]*2)))
}
