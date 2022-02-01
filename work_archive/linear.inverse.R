# linear.inverse <- function(alpha,beta,delta,theta,x){
#   
#   -(delta/(3*theta)) - (1/ (3*2^(1/3)*theta)) * 
#     ((27*alpha*theta^2 - 9*beta*delta*theta+2*delta^3 + 
#     sqrt(4*(3*beta*theta*-delta^2)^3 + (27*alpha*theta^2-9*beta*delta*theta+2*delta^3-27*theta^2*x)^2)-
#     27*theta^2*x)^(1/3)) + 
#     
#     (2^(1/3)*(3*beta*theta-delta^2)) / 
#     (3*theta * (27*alpha*theta^2 - 9*beta*delta*theta+2*delta^3 + sqrt(4*(3*beta*theta-delta^2)^3 +
#                                                                          (27*alpha*theta^2-9*beta*theta+
#                                                                             2*delta^3-27*theta^2*x)^2)-
#                   27*theta^2*x)^(1/3))
# 
#   
# }

linear.inverse <- function(alpha,beta,x){
  ((x - beta) / (alpha) )^(1/3)
  
}

linear.forward <- function(alpha,beta,x){
  alpha*x^3 + beta
  
}
