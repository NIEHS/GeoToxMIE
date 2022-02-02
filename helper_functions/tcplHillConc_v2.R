#-------------------------------------------------------------------------------
# tcplHillConc_v2: Calculate the concentration for a given value, gives 
# concentration in regular space
#-------------------------------------------------------------------------------


tcplHillConc_v2 <- function(val, tp, ga, gw) {
  
  ga * (tp/val -1)^(-1/gw)
  
}

#-------------------------------------------------------------------------------