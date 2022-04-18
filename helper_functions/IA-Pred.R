IA.pred <- function(Ci,tp,AC50,Emax,slope = 1){
  # Independent action prediction from the tcpl Hill function
  
  E <- tcplHillVal_v2(Ci,tp,AC50,slope)
  
  Emax*(1 - prod( (1 - E/Emax) ))
  # 
  
}