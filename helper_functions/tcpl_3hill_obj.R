tcpl.obj.free.slope <- function (p, conc, resp) 
{

    mu <- p[1]/(1 + 10^( p[3] * (p[2] - conc) ))
  

    n <- 4 
    err <- exp(p[n])
    return(sum(dt((resp - mu)/err, df = n, log = TRUE) -
                 log(err)))
}
