# Plotting the concentration-response curves of our data
library(ggplot2)
library(scales)
library(reshape2)

# The concentration-response uses log10 concentration
# We use the dataframe in the function to include the corresponding regular 
# space X
# We plot the x-axis using a log10 scaling 

tp.mean <- get(load("/Volumes/SHAG/GeoTox/data/tp_mean_example.RData"))
AC50.mean <- get(load("/Volumes/SHAG/GeoTox/data/AC50_mean_example.RData"))
slope.mean <- get(load("/Volumes/SHAG/GeoTox/data/slope_mean_example.RData"))

# If we want to include SD somehow in this example
tp.sd <- get(load("/Volumes/SHAG/GeoTox/data/tp_sd_example.RData"))
AC50.sd <- get(load("/Volumes/SHAG/GeoTox/data/AC50_sd_example.RData"))
slope.sd <- get(load("/Volumes/SHAG/GeoTox/data/slope_sd_example.RData"))

myfun <- function(x){
  val <- tcplHillVal(logX,tp.mean[x],AC50.mean[x],slope.mean[x])
  # Testing out the inverse function
  inv.val <- tcplHillConc(val,tp.mean[x],AC50.mean[x],slope.mean[x])
  
  d <- data.frame("logX"= logX, "x"=X,"y"=val,"inv.y" = inv.val)
  return(d)
}


# Same sequence - one is log10 and other is regular space
logX <- log10(seq(10^-3,10^3,length.out = 10000))
X <-seq(10^-3,10^3,length.out = 10000)


val <- lapply(1:41,myfun)

# This melts the list into a dataframe 
m <- melt(val,measure.vars = "y")

ggplot(m,aes(x,value,color = as.factor(L1)))+geom_line()+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
