library(Rcpp)
library(RcppEigen)
sourceCpp("calculate_v_and_a.cpp")
sourceCpp("calculate_HBIC.cpp")
sourceCpp("calculate_initial.cpp")
sourceCpp("calculate_ck.cpp")


LLA = function(X,Y, LL = 100, h = 1, kappa_u = 3/4){
  
  #########################initialize
  n = nrow(X)
  p = ncol(X)
  
  beta_tilde = rep(0,p)#c(rep(sqrt(3),3),rep(0,p-3))
  
  #Rcpp
  ans = calculate_initial(X,Y)
  # ans = rep(0,p)
  # for(i in 1:n){
  #   for(j in 1:n){
  #     if(i != j){
  #       ans = ans + L_prime(Y[i]-Y[j])*(X[i,]-X[j,])
  #     }
  #   }
  # }
  # ans = ans/(n*(n-1))
  
  lambda = (0.1^(1/(LL-1)))^(0:(LL-1))*max(ans) 
  lambda = lambda[floor(0.4*LL):floor(0.7*LL)] ###
  LL = length(lambda)
  solutions = array(0, c(4,LL,p))
  
  #Rcpp
  c_k = calculate_ck(X)
  # c_k = rep(0,p)
  #   for(i in 1:n){
  #     for(j in 1:n){
  #       if(i != j){
  #         c_k = c_k + (X[i,]-X[j,])^2
  #       }
  #     }
  #   }
  # c_k = c_k/(n*(n-1))

  ####################################################
  for( s in 2:4){
    for(l in 1:LL){
      w_k = map_dbl(solutions[s-1,l,1:p],penalty_prime,lambda[l])
      
      beta_bar = rep(0,p)
      d = rep(0,p)
      
      iter = 1
      while(T){
        #Rcpp
        a = calculate_v_and_a(X,Y,beta_tilde)
        # v = matrix(0,n,n)
        # for(i in 1:n){
        #   for(j in 1:n){
        #     v[i,j] = Y[i]-Y[j]-c((X[i,]-X[j,])%*%beta_tilde)
        #   }
        # }
        # 
        # a = rep(0,p)
        # for(i in 1:n){
        #   for(j in 1:n){
        #     if(i != j){
        #       a = a + L_prime(v[i,j])*(X[i,]-X[j,])
        #     }
        #   }
        # }
        # a = -a/(n*(n-1))
        
        ans1 = beta_tilde-h*a/(2*c_k*kappa_u)
        ans2 = h*w_k/(2*c_k*kappa_u)
        ans3 = as.numeric(abs(ans1)>ans2)
        beta_bar = sign(ans1)*(abs(ans1)-ans2)*ans3
        d = beta_bar - beta_tilde
        beta_tilde = beta_bar
        
        iter = iter + 1
        if(max(d^2)<1e-6 | iter>100){break}
      }
      solutions[s,l,] =  beta_tilde
      #cat('s: ',s,';l: ',l,' ','\r')
    }
  }
  
  #Rcpp
  HBIC = calculate_HBIC(solutions[4,,],X,Y)
  # HBIC = apply(solutions[4,,],1,function(b,X,Y){
  #   ans = 0
  #   for(i in 1:n){
  #     for(j in 1:n){
  #       if(i != j){
  #         ans = ans + L_h(Y[i]-Y[j]-c((X[i,]-X[j,])%*%b))
  #       }
  #       
  #     }
  #   }
  #   ans = ans/(n*(n-1))
  #   return(log(ans)+log(log(n))*log(p)/n*sum(b!=0))
  # },X,Y)
  return(solutions[4,which.min(HBIC),])
}

# s1 = Sys.time()
# beta_hat = LLA(X,Y)
# s2 = Sys.time()
# s2-s1
# beta_hat