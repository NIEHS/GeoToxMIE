GCA.obj <- function(effic,Ci,tp,AC50){
  # Generalized concentration addition objective function
  # Trying to find the optimal efficacy value, E
  # Based on a regular space AC50 and concentrations 
  
  # Solving for the efficacy on the natural log-scale. This allows for 
  # better precision in the low values, e.g. 1 x 10-5
  E <- exp(effic)
  ECi <- tcplHillConc_v2(E,tp,AC50,rep(1,length(tp)))
  gca.val <- sum(Ci/ECi,na.rm = F)
  val <- (gca.val - 1)^2
  # print(val)
  return(val)
  
}
