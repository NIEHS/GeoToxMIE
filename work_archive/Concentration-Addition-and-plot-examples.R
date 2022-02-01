# Plotting the concentration-response curves of our data
library(ggplot2)
library(scales)
library(reshape2)


source("/Volumes/messierkp/Projects/AEP-AOP/GeoToxMIE/tcplHillVal.R", echo=FALSE)
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


# Same sequence - one is log10 and other is regular space
logX <- log10(seq(10^-3,10^3,length.out = 10000))
X <-seq(10^-3,10^3,length.out = 10000)



# Equal contributions
Pi <- matrix(1/41,nrow = 41,ncol = 10^4)
ECi <- matrix(NA,nrow = 41,ncol = 10^4)
for (i in 1:41){
  ECi[i,] <- tcplHillVal(logX,tp.mean[i],AC50.mean[i],slope.mean[i])
  
}

Equal.ECmix <-   (colSums(Pi/ECi))^-1

d <- data.frame("logX"= logX, "x"=X,"ECmix"=Equal.ECmix)

ggplot(d,aes(x,ECmix))+geom_line()+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))


for (i in 1:41){
  d <- cbind(d, ECi[i,])
  colnames(d)[i+3] <- paste0("Chem",i) 
  
}

m <- melt(d, measure.vars = 3:44)

ggplot()+geom_line(data = m[m$variable!="ECmix",],aes(x,value,color = variable))+
  geom_line(data = m[m$variable=="ECmix",],aes(x,value),size = 1.5)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme(legend.position = "none")+ggtitle("Equal Proportions- BACKHAUS-ALTENBURGER model") 



### Varying proportions

# some pre-determined proportions - 5 chemicals, all others are 0
p1 <- c(0.5,0.25,0.1,0.1,0.05)
# random chemicals that are present
set.seed(1)
p1.loc <- sample(c(rep(1,5),rep(0,36)),41) %>% as.logical()
Pi <- rep(0,41)
Pi[p1.loc] <- p1
Pi <- replicate(10^4,Pi)
ECi <- matrix(NA,nrow = 41,ncol = 10^4)
for (i in 1:41){
  ECi[i,] <- tcplHillVal(logX,tp.mean[i],AC50.mean[i],slope.mean[i])
  
}

Vary1.ECmix <-   (colSums(Pi/ECi))^-1


d <- data.frame("logX"= logX, "x"=X,"ECmix"=Vary1.ECmix)

for (i in 1:41){
  d <- cbind(d, ECi[i,])
  colnames(d)[i+3] <- paste0("Chem",i) 
  
}

m <- melt(d, measure.vars = 3:44)

ggplot()+geom_line(data = m[m$variable!="ECmix",],aes(x,value,color = variable))+
  geom_line(data = m[m$variable=="ECmix",],aes(x,value),size = 1.5)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme(legend.position = "none")+ggtitle("Varying Proportions- 5 chems- v01- BACKHAUS-ALTENBURGER model") 



# random chemicals that are present
set.seed(100)
p1.loc <- sample(c(rep(1,5),rep(0,36)),41) %>% as.logical()
Pi <- rep(0,41)
Pi[p1.loc] <- p1
Pi <- replicate(10^4,Pi)
ECi <- matrix(NA,nrow = 41,ncol = 10^4)
for (i in 1:41){
  ECi[i,] <- tcplHillVal(logX,tp.mean[i],AC50.mean[i],slope.mean[i])
  
}

Vary1.ECmix <-   (colSums(Pi/ECi))^-1


d <- data.frame("logX"= logX, "x"=X,"ECmix"=Vary1.ECmix)

for (i in 1:41){
  d <- cbind(d, ECi[i,])
  colnames(d)[i+3] <- paste0("Chem",i) 
  
}

m <- melt(d, measure.vars = 3:44)

ggplot()+geom_line(data = m[m$variable!="ECmix",],aes(x,value,color = variable))+
  geom_line(data = m[m$variable=="ECmix",],aes(x,value),size = 1.5)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  theme(legend.position = "none")+ggtitle("Varying Proportions- 5 chems- v02- BACKHAUS-ALTENBURGER model") 

