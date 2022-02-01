# Plotting the concentration-response curves of our data
library(ggplot2)
library(scales)
library(reshape2)
library(gam)
##### source some help functions #####

## tcpl 3 parameter hill function
tcplHillVal <- function(logc, tp, ga, gw, bt = 0) {
  
  bt + (tp - bt)/(1 + 10^((ga - logc)*gw))
  
}

## tcpl 3 parameter hill function inverse
tcplHillConc <- function(val, tp, ga, gw, bt = 0) {
  
  ga - log10((tp - bt)/(val - bt) - 1)/gw
  
}

tcplHillVal2 <- function(c,tp,ga,gw,bt = 0){
  
  bt + (tp - bt) / (1 + (ga / c)^gw)
  
}


tcplHillConc2 <- function(e,tp,ga,gw){
  
  ga * (tp/e -1)^(-1/gw)
  
}
# We
Hill_two <- function(x,alpha,beta){
beta * x / (alpha + x)
}

Hill_two_inv <- function(E,alpha,beta){
  beta*E / (alpha - E)
}

linear.inverse <- function(x,alpha,beta){
  ((x - beta) / (alpha) )^(1/3)
  
}

linear.forward <- function(x,alpha,beta){
  alpha*x^3 + beta
  
}


GCA.obj <- function(E,Ci,tp,AC50,slope){
  
  ECi <- tcplHillConc2(E,tp,AC50,slope)
  idx.na <- is.na(ECi)
  ECi[idx.na] <- tp[idx.na]
  gca.val <- sum(Ci/ECi,na.rm = F)
  # print(sum(is.na(ECi)))
  val <- (gca.val - 1)^2
  # print(c(gca.val,val))
  
  return(val)
  
}

inv.tanh <- function(E,alpha,beta){
  atanh(E/alpha)/beta
}


adjust.obj <- function(param.new,param.old,x,w = rep(1,length(x))){
  # TRY THE SAME THING BUT WITH INVERSE FUNCTIONS INSTEAD
  # param.new = c(tp,AC50)
  # param.old = c(tp,AC50,slope)
  val.old <- tcplHillVal2(x,param.old[1],param.old[2],param.old[3])
  val.new <- tcplHillVal2(x,param.new[1],param.new[2],1)
  
  # val.new <- tcplHillVal2(x,tp,param.new[1],param.new[2])
   # val.new <- Hill_two(x,param.new[1],param.new[2])
  # val.new <- param.new[1] + param.new[2]*x + param.new[3]*x^2 + param.new[4]*x^3
  # val.new <- param.new[1]*tanh(param.new[2]*x)
  # val.new <- linear.forward(x,param.new[1],param.new[2])
  sum((w*(val.old - val.new))^2)
}
### Optimization appraoch
# Goal: Solve for E 
# 1) Initialize a value of E
# 2) Given the concentrations of compounds i = 1,..., n
# 3) Use inverse function to calculate ECxi 
# 4) Calculate the LHS of CA equation ( inv(sum(pi/ECxi)) = 1)
# 5) iterate until converges


tp.mean <- get(load("/Volumes/SHAG/GeoTox/data/archived_outputs/tp_mean_example.RData"))
AC50.mean <- get(load("/Volumes/SHAG/GeoTox/data/AC50_mean_example.RData"))
slope.mean <- get(load("/Volumes/SHAG/GeoTox/data/slope_mean_example.RData"))

X <- 10^seq(log10(10^-3),log10(10^3),length.out = 10^3)

# adjust.obj <- function(param.new,param.old,x){
  # param.new = c(tp,AC50)
  # param.old = c(tp,AC50,slope)

gca.max <- tcplHillVal2(X[length(X)],tp.mean,10^AC50.mean,slope.mean) %>% max()
GCA.val <- rep(NA,length(X))
params.new <- matrix(NA,length(tp.mean),2)
obj.val <- rep(NA,length(tp.mean))
for (i in 1:length(tp.mean)){
  print(c(i))
  # Ci= sort(rnorm(41,mean = logX[i]))

  param.init <- c(tp.mean[i],10^AC50.mean[i])
  # param.init <- rep(0,2)
  # tp.new <- gca.max
  params.slope.one <- optim(par = param.init,fn = adjust.obj,
                  param.old = c(tp.mean[i],10^AC50.mean[i],slope.mean[i]),
                  x = X,
                  w = X/sum(X),
                  method = "Nelder-Mead",
                  control = list(trace=1))
  params.new[i,] <- params.slope.one$par
  # x.design <- cbind(rep(1,length(X)),X,X^2,X^3)
  # y.pred <- x.design %*% params.new[i,]
  # mdl <- lm(y.pred ~ x.design[,2:4])
  # mdl.coef <- coef(mdl)
  obj.val[i] <- params.slope.one$value
}


# df1 <- data.frame("E" = E, "ECmix" = ECmix)
# df2 <- data.frame("E" = E.pred,"C" = conc.total.i$ECmix)
# ggplot()+geom_line(data = df1,aes(ECmix,E))+
#   geom_point(data = df1,aes(ECmix,E))+
#   geom_point(data = df2,aes(C,E),color = "red")

for (i in 1:41){
  mdl1A <- tcplHillVal2(X,tp.mean[i],10^AC50.mean[i],slope.mean[i])
  mdl1B <- tcplHillVal2(X,params.new[i,1],params.new[i,2],1)
  # mdl1B <- Hill_two(X,params.new[1,1],params.new[1,2])
  # mdl1B <- params.new[i,1] + params.new[i,2]*X + params.new[i,3]*X^2 +
  #   params.new[i,4]*X^3
  # mdl1B <- params.new[i,1]*tanh(params.new[i,2]*X)
  # mdl1B <- linear.forward(X,params.new[i,1],params.new[i,2])
  df1 <- data.frame("X" = X, "Y1"= mdl1A, "Y2" = mdl1B)
  p <- ggplot(df1)+geom_line(aes(X,Y1),color = "red")+geom_line(aes(X,Y2),color = "blue")+
  scale_x_log10() 
  print(p)
}

