# Plotting the concentration-response curves of our data
library(ggplot2)
library(scales)
library(reshape2)
library(mixtox)
##### source some help functions #####

## tcpl 3 parameter hill function
tcplHillVal <- function(logc, tp, ga, gw, bt = 0) {
  
  bt + (tp - bt)/(1 + 10^((ga - logc)*gw))
  
}

## tcpl 3 parameter hill function inverse
tcplHillConc <- function(val, tp, ga, gw, bt = 0) {
  
  ga - log10((tp - bt)/(val - bt) - 1)/gw
  
}

# Weibull
# E = 1 − exp(− exp(α + β log(c)))
Weibull <- function(logc,alpha,beta){
  
  1 - exp(-exp(alpha + beta*logc))
}

# Logit
# E = (1 + exp(−α − β log(c)))−1

Logit <- function(logc,alpha,beta){
  
  (1 + exp(-alpha - beta *logc))^-1
}

# Inverse Weibull
# c = 10(ln(− ln(1−E))−α)/β
inv.Weibull <- function(E,alpha,beta){
  ((log(-log(1-E))-alpha)/beta)
}

# Inverse Logit
# 10^(ln(E/(1−E))−α)/β
inv.logit <- function(E,alpha,beta){
  (log(E/(1-E))-alpha)/beta
}

# Concentration addition objective function to find mixture effect


CA.obj <- function(ECmix,E,Ci,alpha,beta,model){
  
  pi <- Ci / sum(Ci)
  ECi <- rep(NA,length(model))
  for (i in 1:length(model)){
    if (model[i]=="Logit"){
      ECi[i] <- 10^inv.logit(E,alpha[i],beta[i])
    }else if(model[i] == "Weibull"){
     ECi[i] <- 10^inv.Weibull(E,alpha[i],beta[i])  
    }
  }
  ca.val <- sum((pi * ECmix)/ECi,na.rm = F)
  # print(sum(is.na(ECi)))
  val <- (ca.val - 1)^2
  
   # print(c(ECmix,val))
  return(val)
  

  
}

GCA.obj <- function(E,Ci,alpha,beta,model){
  
  ECi <- rep(NA,length(model))
  for (i in 1:length(model)){
    if (model[i]=="Logit"){
      ECi[i] <- 10^inv.logit(E,alpha[i],beta[i])
    }else if(model[i] == "Weibull"){
      ECi[i] <- 10^inv.Weibull(E,alpha[i],beta[i])  
    }
  }
  ca.val <- sum(Ci/ECi,na.rm = F)
  # print(sum(is.na(ECi)))
  val <- (ca.val - 1)^2
  
  # print(c(ECmix,val))
  return(val)
  
  
  
}



##### END help functions #####

# Mixtox example
model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
mixtox.pred <- caPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))

df <- data.frame("ca5" = mixtox.pred$ca[1,],"ca50" = mixtox.pred$ca[2,],
                 "e" = mixtox.pred$e)

chem.names <- antibiotox$sgl$param %>% rownames()
xi <- matrix(NA,nrow = 7,ncol = 12)
for (i in 1:7){
  xi[i,] <- antibiotox[[i]]$x
}

x <- colSums(xi)
logc <- log10(x)
E <- seq(1e-4,1-1e-4,length.out = 100)
E.val <- rep(NA,length(logc))
for (i in 1:length(logc)){
  print(c(i))
  Ci= (c(antibiotox$PAR$x[i],
              antibiotox$SPE$x[i],
              antibiotox$KAN$x[i],
              antibiotox$STR$x[i],
              antibiotox$DIH$x[i],
              antibiotox$GEN$x[i],
              antibiotox$NEO$x[i]))
  
  ECmix <-  rep(NA,length(E))
  for (j in 1:length(E)){
    mixture.result <- optimize(f = CA.obj, interval = c(-20,20),
                               # Ci= rep(X[i],41),
                               E = E[j],
                               Ci = Ci,
                               alpha = param[,1],
                               beta = param[,2],
                               model = antibiotox$sgl$model)
    ECmix[j] <- mixture.result$minimum 
  }
  

  mdl <- gam::gam(E ~ s(ECmix,df = 98))
  
  pi <- Ci/sum(Ci)
  conc.total.i <- data.frame("ECmix"= sum(Ci * pi))
  E.val[i] <- predict(mdl,newdata =  conc.total.i)
  

}

df <- data.frame("X" = x,"Y" = E.val)

ggplot()+geom_line(data = df,aes(X,Y),size =1.5)+geom_point(data = df,aes(X,Y))+
  geom_point(data = data.frame("X"=(antibiotox$PAR$x),"Y" = antibiotox$PAR$y[,1]),
            aes(X,Y),color = "Red")+
  geom_line(data = data.frame("X"=(antibiotox$PAR$x),"Y" = antibiotox$PAR$y[,1]),
             aes(X,Y),color = "Red")+
  geom_point(data = data.frame("X"=(antibiotox$SPE$x),"Y" = antibiotox$SPE$y[,1]),
             aes(X,Y),color = "Green")+
  geom_line(data = data.frame("X"=(antibiotox$SPE$x),"Y" = antibiotox$SPE$y[,1]),
             aes(X,Y),color = "Green")+
  geom_point(data = data.frame("X"=(antibiotox$KAN$x),"Y" = antibiotox$KAN$y[,1]),
             aes(X,Y),color = "coral3")+
  geom_line(data = data.frame("X"=(antibiotox$KAN$x),"Y" = antibiotox$KAN$y[,1]),
             aes(X,Y),color = "coral3")+
  geom_point(data = data.frame("X"=(antibiotox$STR$x),"Y" = antibiotox$STR$y[,1]),
             aes(X,Y),color = "tomato")+
  geom_line(data = data.frame("X"=(antibiotox$STR$x),"Y" = antibiotox$STR$y[,1]),
             aes(X,Y),color = "tomato")+
  geom_point(data = data.frame("X"=(antibiotox$DIH$x),"Y" = antibiotox$DIH$y[,1]),
             aes(X,Y),color = "bisque1")+
  geom_line(data = data.frame("X"=(antibiotox$DIH$x),"Y" = antibiotox$DIH$y[,1]),
             aes(X,Y),color = "bisque1")+
  geom_point(data = data.frame("X"=(antibiotox$GEN$x),"Y" = antibiotox$GEN$y[,1]),
             aes(X,Y),color = "plum3")+
  geom_line(data = data.frame("X"=(antibiotox$GEN$x),"Y" = antibiotox$GEN$y[,1]),
             aes(X,Y),color = "plum3")+
  geom_point(data = data.frame("X"=(antibiotox$NEO$x),"Y" = antibiotox$NEO$y[,1]),
             aes(X,Y),color = "slateblue1")+
  geom_line(data = data.frame("X"=(antibiotox$NEO$x),"Y" = antibiotox$NEO$y[,1]),
             aes(X,Y),color = "slateblue1")+
  scale_x_log10()


### GCA Testing

GCA.val <- rep(NA,length(logc))
for (i in 1:length(logc)){
  print(c(i))
  Ci= (c(antibiotox$PAR$x[i],
         antibiotox$SPE$x[i],
         antibiotox$KAN$x[i],
         antibiotox$STR$x[i],
         antibiotox$DIH$x[i],
         antibiotox$GEN$x[i],
         antibiotox$NEO$x[i]))
  
  ECmix <-  rep(NA,length(E))
    mixture.result <- optimize(f = GCA.obj, interval = c(10^-20,1),
                               # Ci= rep(X[i],41),
                               Ci = Ci,
                               alpha = param[,1],
                               beta = param[,2],
                               model = antibiotox$sgl$model)


  GCA.val[i] <- mixture.result$minimum

  
}

df <- data.frame("X" = x,"CA" = E.val,"GCA" = GCA.val)


ggplot()+geom_line(data = df,aes(X,CA),size =1.5)+geom_point(data = df,aes(X,CA))+
  geom_line(data = df,aes(X,GCA),size =1.5,color = "goldenrod")+geom_point(data = df,aes(X,GCA),color = "goldenrod")+
  geom_point(data = data.frame("X"=(antibiotox$PAR$x),"Y" = antibiotox$PAR$y[,1]),
             aes(X,Y),color = "Red")+
  geom_line(data = data.frame("X"=(antibiotox$PAR$x),"Y" = antibiotox$PAR$y[,1]),
            aes(X,Y),color = "Red")+
  geom_point(data = data.frame("X"=(antibiotox$SPE$x),"Y" = antibiotox$SPE$y[,1]),
             aes(X,Y),color = "Green")+
  geom_line(data = data.frame("X"=(antibiotox$SPE$x),"Y" = antibiotox$SPE$y[,1]),
            aes(X,Y),color = "Green")+
  geom_point(data = data.frame("X"=(antibiotox$KAN$x),"Y" = antibiotox$KAN$y[,1]),
             aes(X,Y),color = "coral3")+
  geom_line(data = data.frame("X"=(antibiotox$KAN$x),"Y" = antibiotox$KAN$y[,1]),
            aes(X,Y),color = "coral3")+
  geom_point(data = data.frame("X"=(antibiotox$STR$x),"Y" = antibiotox$STR$y[,1]),
             aes(X,Y),color = "tomato")+
  geom_line(data = data.frame("X"=(antibiotox$STR$x),"Y" = antibiotox$STR$y[,1]),
            aes(X,Y),color = "tomato")+
  geom_point(data = data.frame("X"=(antibiotox$DIH$x),"Y" = antibiotox$DIH$y[,1]),
             aes(X,Y),color = "bisque1")+
  geom_line(data = data.frame("X"=(antibiotox$DIH$x),"Y" = antibiotox$DIH$y[,1]),
            aes(X,Y),color = "bisque1")+
  geom_point(data = data.frame("X"=(antibiotox$GEN$x),"Y" = antibiotox$GEN$y[,1]),
             aes(X,Y),color = "plum3")+
  geom_line(data = data.frame("X"=(antibiotox$GEN$x),"Y" = antibiotox$GEN$y[,1]),
            aes(X,Y),color = "plum3")+
  geom_point(data = data.frame("X"=(antibiotox$NEO$x),"Y" = antibiotox$NEO$y[,1]),
             aes(X,Y),color = "slateblue1")+
  geom_line(data = data.frame("X"=(antibiotox$NEO$x),"Y" = antibiotox$NEO$y[,1]),
            aes(X,Y),color = "slateblue1")+
  scale_x_log10()+
  ggtitle("Testing CA and GCA Mixture Effect Predictions")+
  xlab("Concentration")+
  ylab("Effect")
