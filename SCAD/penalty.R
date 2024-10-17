penalty = function(t,lambda,a=3.7){
  if(abs(t)<lambda){
    return(lambda*abs(t))
  }else if(abs(t)<=a*lambda){
    return(
      1/(a-1)*(a*lambda*abs(t)-(t^2+lambda^2)/2)
    )}
  else{
      return(
        (a+1)*lambda^2/2
      )
  }
}

penalty_prime = function(t,lambda,a=3.7){
  if(0<=t & t<lambda){
    return(lambda)
  }else if(0>t & t>-lambda){
    return(-lambda)
  }else if(lambda<=t & t<=a*lambda){
    return(
      1/(a-1)*(a*lambda - t)
    )
  }else if(-lambda>=t & t>=-a*lambda){
    return(
      1/(a-1)*(-a*lambda - t)
    )
  }else{
    return(0)
  }
}

