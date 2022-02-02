IA.pred <- function(Ci,tp,AC50,slope = 1){
  # Independent action prediction from the tcpl Hill function
  
  E <- tcplHillVal_v2(Ci,tp,AC50,slope)
  
  1 - prod( (1 - E) )
  
}
