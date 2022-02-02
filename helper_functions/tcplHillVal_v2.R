#-------------------------------------------------------------------------------
# tcplHillVal_v2: Calculate the efficacy for a given concentration (regular space concentration)
# Regular space AC50 too
#-------------------------------------------------------------------------------


tcplHillVal_v2 <- function(c, tp, ga, gw, bt = 0) {
  
  bt + (tp - bt) / (1 + (ga / c)^gw)
  
}
