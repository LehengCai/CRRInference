library(MASS)

Gen_data = function(n, p, beta, X_type = 1, error_type = 1){
  if(X_type == 1){
    rho = 0.5
    temp = rho^seq(0,p-1,1)
    mu = rep(0,p)
    Sigma = toeplitz(temp)
    X = mvrnorm(n, mu, Sigma) 
  }else if(X_type == 2){
    rho = 0.48
    Sigma = matrix(0,p,p)
    diag(Sigma) = rep(1,p)

    Sigma[row(Sigma) == (col(Sigma) - 1)] <- rho
    Sigma[row(Sigma) == (col(Sigma) + 1)] <- rho
    mu = rep(0,p)
    X = mvrnorm(n, mu, Sigma) 
  }else if(X_type == 3){
    Sigma = diag(rep(1,p))
    mu = rep(0,p)
    X = mvrnorm(n, mu, Sigma)
  }
  
  if(error_type == 1){
    epsilon = rnorm(n,0,1)
  }else if(error_type == 2){
    epsilon = rt(n,df = 4)/sqrt(2)
  }else if(error_type == 3){
    epsilon = (rgamma(n,4,1)-2)/2
  }else if(error_type == 4){
    epsilon = rt(n,df = 3)
  }else if(error_type == 5){
    epsilon = rcauchy(n)
  }else if(error_type == 6){
    z = rbinom(n,1,prob=0.95)
    epsilon = rnorm(n,0,1)*z + rnorm(n,0,100)*(1-z)
  }
  
  Y = X%*%beta+epsilon
  Y = as.vector(Y)
  return(list(X=X,Y=Y))
}

# n = 100
# p = 50
# beta = c(rep(sqrt(3),3),rep(0,p-3))
# Data = Gen_data(n,p,beta)
# X = Data$X
# Y = Data$Y
