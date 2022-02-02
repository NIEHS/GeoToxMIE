#-------------------------------------------------------------------------------
# tcplHillVal: Calculate the efficacy for a given concentration (log-concentration)
#-------------------------------------------------------------------------------


tcplHillVal <- function(logc, tp, ga, gw, bt = 0) {
  
  bt + (tp - bt)/(1 + 10^((ga - logc)*gw))
  
}
