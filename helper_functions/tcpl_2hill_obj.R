tcpl.obj.fixed.slope <- function (p, conc, resp) 
{

  mu <- p[1]/(1 + 10^( (p[2] - conc) ))
  

    n <- 3 
    err <- exp(p[n])
    return(sum(dt((resp - mu)/err, df = n, log = TRUE) -
                 log(err)))
}
