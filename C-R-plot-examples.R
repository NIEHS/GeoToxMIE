# Plotting the concentration-response curves of our data
library(ggplot2)
library(scales)
library(reshape2)

# The concentration-response uses log10 concentration
# We use the dataframe in the function to include the corresponding regular 
# space X
# We plot the x-axis using a log10 scaling 

tp.mean <- get(load("/Volumes/SHAG/GeoTox/data/tp_mean_example.RData"))

function(x){
  val <- tcplHillVal(logX,tp.mean[x],AC50.mean[x],slope.mean[x])
  d <- data.frame("x"=X,"y"=val)
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
