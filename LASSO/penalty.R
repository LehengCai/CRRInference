penalty = function(t,lambda){
    return(lambda*abs(t))
}

penalty_prime = function(t,lambda){
  if(t>=0){
    return(lambda)
  }else{
    return(-lambda)
  }
}

