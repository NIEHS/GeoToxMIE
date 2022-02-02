GCA.obj <- function(E,Ci,tp,AC50){
  # Generalized concentration addition objective function
  # Trying to find the optimal efficacy value, E
  # Based on a regular space AC50 and concentrations 
  
  ECi <- tcplHillConc_v2(E,tp,AC50,1)
  gca.val <- sum(Ci/ECi,na.rm = F)
  val <- (gca.val - 1)^2

  return(val)
  
}