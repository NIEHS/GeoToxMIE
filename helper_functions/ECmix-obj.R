ECx.obj <- function(ECmix,E,Ci,tp,AC50){
  # Find the effective concentration of a mixture via an objective function
  # given the concentrations and inverse 
  # A Generalized concentration addition style solution
  # Based on a regular space AC50 and concentrations 
  
  ECi <- tcplHillConc_v2(E,tp,AC50,rep(1,length(tp)))
  
  ECx.val <- sum( pi * ECmix / ECi,na.rm = F)
  val <- (ECx.val - 1)^2
  return(val)
  
}
