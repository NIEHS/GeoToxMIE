my_2hill_fit <- function(conc, resp, log = FALSE, verbose = FALSE) 
{
  
  if (log == TRUE){ # Resposne data is already log-transformed
    logc = conc
  }else{
    logc = log10(conc)
  }

  rmds <- tapply(resp, logc, median)
  
  mmed = rmds[which.max(abs(rmds))]
  mmed_conc <- as.numeric(names(mmed))
  resp_max <- max(resp)
  resp_min <- min(resp)
  logc_min <- min(logc)
  logc_max <- max(logc)
  
  er_est <- if ((rmad <- mad(resp)) > 0) 
    log(rmad)
  else log(1e-32)
  
  # lb <- c(0,logc_min - 2,0.3,-Inf)
  # ub <- c(1.2*resp_max,logc_max +0.5,8,Inf)
  # 
  # param <- c(mmed, mmed_conc - 0.5,1.2,er_est)
  lb <- c(0,logc_min - 2,-Inf)
  ub <- c(1.2*resp_max,logc_max +0.5,Inf)
  
  param <- c(mmed, mmed_conc - 0.5,er_est)
  
  fit <- optim(param, fn = tcpl.obj.fixed.slope,
               method = "L-BFGS-B",
               conc = logc, 
               resp = resp,
               lower = lb,
               upper = ub,
               hessian = TRUE,
               control = list(fnscale = -1,
                              maxit = 10000))
  
  aic <- 2 * length(fit$par) - 2 * fit$value
  
  sds <- sqrt(diag(solve(-fit$hessian)))
  param.names <- c("tp","logAC50","t-error")
  mdl <- list()
  mdl$par <- fit$par
  mdl$sds <- sds
  names(mdl$par) <- param.names
  names(mdl$sds) <- param.names
  mdl$val <- fit$value
  mdl$convergence <- fit$convergence
  mdl$AIC <- aic
  mdl$logc_max <- logc_max
  mdl$logc_min <- logc_min
  mdl$resp_max <- resp_max
  mdl$resp_min <- resp_min

  return(mdl)
  
}
