CA.obj <- function(E,Ci,tp,AC50,slope){
# Concentration addition objective function
# Input:
# E: scalar, the mixture Effect (i.e. response)
# Ci: vector, the concentrations of the chemicals

### Version 1 - Not working
# The mixture effect is common to each chemical
# Mixture.Effect <- rep(E,length(tp))
# # Inverse of 3-parameter hill model to calculate the effect concentration  
# 
# # ECmix
ECxi <- tcplHillConc(Mixture.Effect,tp,AC50,slope)
idx.na <- !is.na(ECxi)
ECxi <- ECxi[idx.na]
p.i <- Ci[idx.na] / sum(Ci[idx.na])
ECmix <- (sum(p.i/ECxi,na.rm = TRUE))^-1

Ci <- p.i * ECmix
# Calculate the left-hand side of the equation
Sigma.val <- Ci/ECxi
Sigma.val[is.infinite(Sigma.val)] <- NA
lhs <- sum(Sigma.val,na.rm = TRUE)
### End Version 1 - Not working

# Version 2 - not working
  
# ECxi <- tcplHillConc(Mixture.Effect,tp,AC50,slope)
# Sigma.val <- Ci/ECxi
# Sigma.val[is.infinite(Sigma.val)] <- NA
# lhs <- sum(Sigma.val,na.rm = TRUE)
# end version 2

  # Version 3
  # pi <- Ci / sum(Ci)
  # ECi <- tcplHillConc(Mixture.Effect,tp,AC50,slope)
  # # idx.na <- !is.na(ECi)
  # ECmix <- sum( (pi/ECi),na.rm = TRUE)^-1
  #  
  # print(ECmix)
  # val <- (ECmix - E)^2 
# Minimize the sum of squared error from 1
val <- (lhs - 1)^2
  
# print(c(val,lhs))
return(val)

}








  