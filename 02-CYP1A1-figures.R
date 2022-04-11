######################################################
# Manuscript Figure S5
# By: Kristin Eccles and Kyle Messier
# Date: April 11th, 2022
# Written in R Version 4.0.2
######################################################

# Fit concentration-response data to the 2-parameter hill function
# Compare to the 3-parameter ICE data
# Check the standard deviation estimates, etc. 

library(scales)
library(reshape2)
library(stringr)
library(dplyr)
library(ggplot2)
library(sjPlot)

### Load in the in-vitro concentration response data from ICE
ice_data <- get(load("/Volumes/SHAG/GeoTox/data/LTEA_HepaRG_CYP1A1_up 41 chems for Kyle 220131.RData"))
unique(ice_data$chnm, ice_data$casn)
hill2.fit <- get(load("/Volumes/SHAG/GeoTox/data/Hill_2param_model_fit.RData"))

# The concentration-response uses log10 concentration
# We use the dataframe in the function to include the corresponding regular 
# space X
# We plot the x-axis using a log10 scaling 
tcplHillVal <- function(logc, tp, ga, gw, bt = 0) {
  
  bt + (tp - bt)/(1 + 10^((ga - logc)*gw))
  
}

# tp.mean <- get(load("/Volumes/SHAG/GeoTox/data/archived_outputs/tp_mean_example.RData"))
# AC50.mean <- get(load("/Volumes/SHAG/GeoTox/data/AC50_mean_example.RData"))
# slope.mean <- get(load("/Volumes/SHAG/GeoTox/data/slope_mean_example.RData"))

# load data
hill2 <- get(load("/Volumes/SHAG/GeoTox/data/Hill_2param_model_fit.RData"))

# Example hill model parameters (from the EPA/ICE database)

tp.mean <- hill2$tp
AC50.mean <- (hill2$logAC50)
AC50.mean.normal <- 10^((hill2$logAC50))

tp.sd<- (hill2$tp.sd)
AC50.sd <- (hill2$logAC50.sd)
AC50.sd.normal <- 10^(hill2$logAC50.sd)

tp.lower<- tp.mean - tp.sd 
AC50.lower <- log10(AC50.mean.normal - AC50.sd.normal)
#AC50.lower[AC50.lower < 0] <- 1

tp.upper<- tp.mean + tp.sd
AC50.upper <- log10(AC50.mean.normal + AC50.sd.normal)

slope <- c(1,1,1,1,1,1,1,1,1,1,
           1,1,1,1,1,1,1,1,1,1,
           1,1,1,1,1,1,1,1,1,1,1,
           1,1,1,1,1,1,1,1,1,1)

myfun_mean <- function(x){
  val <- tcplHillVal(logX,tp.mean[x],AC50.mean[x],slope[x])

  d <- data.frame("logX"= logX, "x"=X,"mean"=val)
  return(d)
}

myfun_lower<- function(x){

  val <- tcplHillVal(logX,tp.lower[x],AC50.lower[x],slope[x])
  
  d <- data.frame("logX"= logX, "x"=X,"lower"=val)
  return(d)
}

myfun_upper<- function(x){
  
  val <- tcplHillVal(logX,tp.upper[x],AC50.upper[x],slope[x])
  
  d <- data.frame("logX"= logX, "x"=X,  "upper" = val)
  return(d)
}


# Same sequence - one is log10 and other is regular space
logX <- log10(seq(10^-1,10^3,length.out = 10000))
X <-seq(10^-1,10^3,length.out = 10000)


val_mean <- lapply(1:41,myfun_mean)
val_lower <- lapply(1:41,myfun_lower)
val_upper <- lapply(1:41,myfun_upper)

# This melts the list into a dataframe 
mean_melt <- melt(val_mean, measure.vars = "mean")
lower_melt <- melt(val_lower, measure.vars = "lower")
upper_melt <- melt(val_upper, measure.vars = "upper")

m <- cbind(mean_melt, lower_melt[,4], upper_melt[,4])
#m <- cbind(mean_melt, lower_melt, upper_melt)
colnames(m) <- c("logX",  "x", "var" , "mean", "id", "lower", "upper")
#m <- m[!duplicated(as.list(m))]

hill2$id <- 1:nrow(hill2)
hill2$chnm <- rownames(hill2)

# replace chem name names that are wrong
#C.I. Disperse Black 6 119-90-4 = 3,3'-Dimethoxybenzidine 
#C.I. Azoic Diazo Component 112 = 92-87-5 = Benzidine
hill2$chnm <- str_replace_all(hill2$chnm, "C.I. Disperse Black 6", "3,3'-Dimethoxybenzidine")
hill2$chnm <- str_replace_all(hill2$chnm, "C.I. Azoic Diazo Component", "Benzidine")


m <- left_join(m, hill2[,c("id", "chnm")], by="id", keep=FALSE)
m$lower[m$lower < 0] <- 0
  
SI_figure<- ggplot()+

  geom_line(data = m, aes(x = x, y = lower), color="red")+
  geom_line(data = m, aes(x = x, y = upper), color="red")+
  geom_line(data = m, aes(x = x, y = mean), color="black", size = 1)+
  #geom_point(data=ice_data, aes(x=logc, y=resp))+
  facet_wrap(vars(chnm), scales= "free")+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme_bw()+
  theme(legend.position = "none",
        text = element_text(size = 16))+
  ylab("Response")+
  xlab("Concentration")
SI_figure

save_plot("SI_dr.tif", SI_figure, width= 50, heigh=40, dpi=300)

