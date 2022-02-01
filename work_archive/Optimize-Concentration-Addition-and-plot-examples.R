# Plotting the concentration-response curves of our data
library(ggplot2)
library(scales)
library(reshape2)


source("/Volumes/messierkp/Projects/AEP-AOP/GeoToxMIE/tcplHillVal.R", echo=FALSE)
source("/Volumes/messierkp/Projects/AEP-AOP/GeoToxMIE/tcplHillConc.R", echo=FALSE)
source("/Volumes/messierkp/Projects/AEP-AOP/GeoToxMIE/concentration_addition.R")
# The concentration-response uses log10 concentration
# We use the dataframe in the function to include the corresponding regular 
# space X
# We plot the x-axis using a log10 scaling 

tp.mean <- get(load("/Volumes/SHAG/GeoTox/data/tp_mean_example.RData"))
AC50.mean <- get(load("/Volumes/SHAG/GeoTox/data/AC50_mean_example.RData"))
slope.mean <- get(load("/Volumes/SHAG/GeoTox/data/slope_mean_example.RData"))

# # If we want to include SD somehow in this example
# tp.sd <- get(load("/Volumes/SHAG/GeoTox/data/tp_sd_example.RData"))
# AC50.sd <- get(load("/Volumes/SHAG/GeoTox/data/AC50_sd_example.RData"))
# slope.sd <- get(load("/Volumes/SHAG/GeoTox/data/slope_sd_example.RData"))


# Same sequence - one is log10 and other is regular space
logX <- log10(seq(10^-3,10^3,length.out = 10^3))
X <-seq(10^-3,10^3,length.out = 10^3)


# Goal: Solve for E 
# 1) Initialize a value of E
# 2) Given the concentrations of compounds i = 1,..., n
# 3) Use inverse function to calculate ECxi 
# 4) Calculate the LHS of CA equation
# 5) Calculate differnce from 1 
# 6) iterate until converges

# 1)
Mixture.Effect <- 1
#CA.obj <- function(E,Ci,tp,AC50,slope){

init.upper <- tcplHillVal(logX[10^3],tp.mean,AC50.mean,slope.mean) %>% 
  max() 
# Equal contributions
E.val <- rep(NA,length(logX))
for (i in 1:length(logX)){
  print(c(i))
  mixture.result <- optim(par=Mixture.Effect,fn=CA.obj,
                          Ci= rep(logX[i],41),
                          tp = tp.mean,
                          AC50 = AC50.mean,
                          slope = slope.mean,
                          method = "Nelder-Mead",
                          control=list(trace=0))
  
  mixture.result <- optimize(f = CA.obj, interval <- c(0,init.upper),
                             Ci= rep(logX[i],41),
                             tp = tp.mean,
                             AC50 = AC50.mean,
                             slope = slope.mean)

  # E.val[i] <- mixture.result$par
  E.val[i] <- mixture.result$minimum
  
}


## 
init.upper <- tcplHillVal(logX[10^3],tp.mean,AC50.mean,slope.mean) %>% 
  max() 

init.lower <- tcplHillVal(logX[1],tp.mean,AC50.mean,slope.mean) %>% 
  min() 
mixture.seq <- seq(1e-20,1,length.out = 10^3)

xx.seq <- seq(0,100,length.out = 10^3)
ECmix <- rep(NA, 10^3)
for (i in 1:10^3){
  ECi <- tcplHillConc(mixture.seq[i],tp.mean,AC50.mean,slope.mean)
  # ECi[is.na(ECi)] <- init.lower
  ECmix[i] <- (sum(rep(1/41,41)/ECi,na.rm = TRUE))^-1
  
}

df <- data.frame("X" = 10^ECmix, "logX" = ECmix, "Y" = mixture.seq)

ggplot(df,aes(logX,Y))+ geom_line()+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))

pi <- matrix(1/41,nrow = 10^3,ncol = 41)
ECmix <- (rowSums(pi/ECxi,na.rm = T)^-1)

