setwd('C:/Users/Cai/Desktop/Research/CRRInference/major_revision/Github/LASSO')
library(Rcpp)
library(RcppEigen)
library(stringr)

source('Gen_data.R')
source('kernel.R')
source('penalty.R')
source('estimation_HBIC.R')
source('inference.R')

test = function(xxx,n,p,beta,alpha,h,B,X_type = 1, error_type = 1){
  cat('Running the ',xxx,'-th iteration \r')
  Data = Gen_data(n,p,beta,X_type, error_type)
  X = Data$X
  Y = Data$Y
  beta_hat = LLA(X,Y,h = h)
  result = inference(X,Y,beta_hat,beta,alpha,B=B)
  return(result)
}
##############################################
MCMC = 500

set.seed(123)
B = 500
h = 1
n = 100
p = 50
X_type = 1
error_type = 1
alpha = 0.05
beta = c(rep(sqrt(3),3),rep(0,p-3))

s1 = Sys.time()
f = sapply(1:MCMC,test,n,p,beta,alpha,h,B,X_type,error_type)
s2 = Sys.time()
s2-s1
apply(f,1,mean)

######################################################################################
setwd('C:/Users/Cai/Desktop/Research/CRRInference/major_revision/Github/SCAD')
library(Rcpp)
library(RcppEigen)
library(stringr)

source('Gen_data.R')
source('kernel.R')
source('penalty.R')
source('estimation_HBIC.R')
source('inference.R')

test = function(xxx,n,p,beta,alpha,h,B,X_type = 1, error_type = 1){
  cat('Running the ',xxx,'-th iteration \r')
  Data = Gen_data(n,p,beta,X_type, error_type)
  X = Data$X
  Y = Data$Y
  beta_hat = LLA(X,Y,h = h)
  result = inference(X,Y,beta_hat,beta,alpha,B=B)
  return(result)
}
##############################################
MCMC = 500

set.seed(123)
B = 500
h = 1
n = 100
p = 50
X_type = 1
error_type = 1
alpha = 0.05
beta = c(rep(sqrt(3),3),rep(0,p-3))

s1 = Sys.time()
f = sapply(1:MCMC,test,n,p,beta,alpha,h,B,X_type,error_type)
s2 = Sys.time()
s2-s1
apply(f,1,mean)