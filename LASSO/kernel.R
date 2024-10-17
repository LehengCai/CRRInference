library(purrr)

K = function(u){
  if(u>=-1 & u<=1){
    return( 3/4*(1-u^2) )
  }else{
    return(0)
  }
}
K_h = function(u,h=1){
  return(1/h*K(u/h))
}

L = function(u){
  if((u >= -1) & (u <= 1)){
    return( 3/4*u^2-u^4/8+3/8 )
  }else if (u>1){
    return(u)
  }else{
    return(-u)
  }
}
L_h = function(u,h=1){
  return(h*L(u/h))
}

L_prime = function(u){
  if((u >= -1) & (u <= 1)){
    return( 3/2*u-u^3/2)
  }else if (u>1){
    return(1)
  }else{
    return(-1)
  }
}
