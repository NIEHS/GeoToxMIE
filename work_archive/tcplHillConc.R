#-------------------------------------------------------------------------------
# tcplHillConc: Calculate the concentration for a given value
#-------------------------------------------------------------------------------

#' @rdname hill_utils
#' @export

tcplHillConc <- function(val, tp, ga, gw, bt = 0) {
  
  val <- rep(val,length(tp))
  
  idx <- val > tp
  
  val[idx] <- tp[idx] - 1e-10
  
  ga - log10((tp - bt)/(val - bt) - 1)/gw

}

#-------------------------------------------------------------------------------