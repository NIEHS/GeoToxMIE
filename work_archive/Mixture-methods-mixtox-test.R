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

tcplHillVal2 <- function(c,tp,ga,gw){
  
  tp  / (1 + (ga / c)^gw)
  
}


tcplHillConc2 <- function(e,tp,ga,gw){

  ga * (tp/e -1)^(-1/gw)
  
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

inv.hill.two <- function(E,alpha,beta){
  # βE/ (α − E)
  (beta*E)/(alpha - E)
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
  gca.val <- sum(Ci/ECi,na.rm = F)
  # print(sum(is.na(ECi)))
  val <- (gca.val - 1)^2
  # print(c(gca.val,val))
  
  return(val)
  
}

IA.pred <- function(Ci,alpha,beta){

  # pi <- Ci / sum(Ci)
  # Ci <- Ci * pi
  E <- rep(NA,length(model))
  for (i in 1:length(model)){
    if (model[i]=="Logit"){
      E[i] <- Logit(Ci[i],alpha[i],beta[i])
    }else if(model[i] == "Weibull"){
      E[i] <- Weibull(Ci[i],alpha[i],beta[i])  
    }
  }
  
  1 - prod( (1 - E) )
  
}


##### END help functions #####


# Example model parameters and model
param <- antibiotox$sgl$param
model <- antibiotox$sgl$model


#

xi <- matrix(NA,nrow = 7,ncol = 12)
for (i in 1:7){
  xi[i,] <- antibiotox[[i]]$x
}

pi <- matrix(NA,nrow = 7, ncol = 12)
for (i in 1:12){
  pi[,i] <- xi[,i] / sum(xi[,i]) 
}



# Concentration totals 
x.total <- colSums(xi)

# Effect sequence used for 2-stage CA approach
E <- seq(1e-4,1-1e-4,length.out = 100)


### GCA Testing

CA.val <- GCA.val <- AVG.val <- IA.val <- rep(NA,length(x.total))

for (i in 1:length(x.total)){
  print(c(i))
  Ci= (c(antibiotox$PAR$x[i],
         antibiotox$SPE$x[i],
         antibiotox$KAN$x[i],
         antibiotox$STR$x[i],
         antibiotox$DIH$x[i],
         antibiotox$GEN$x[i],
         antibiotox$NEO$x[i]))

 #### CA, Concentration Addition 
  ECmix <-  rep(NA,length(E))
  for (j in 1:length(E)){
    mixture.result <- optimize(f = CA.obj, interval = c(-20,20),
                               E = E[j],
                               Ci = Ci,
                               alpha = param[,1],
                               beta = param[,2],
                               model = antibiotox$sgl$model)
    ECmix[j] <- mixture.result$minimum 
  }
  
  conc.total.i <- data.frame("ECmix"= sum(Ci))
  
  idx.which <- which.min(( ECmix - conc.total.i$ECmix)^2)
  # mdl <- gam::gam(E ~ s(ECmix,df = 98))
  # 
  # # pi <- Ci/sum(Ci)
  # CA.val[i] <- predict(mdl,newdata =  conc.total.i)
  CA.val[i] <- E[idx.which]
  
  ### End CA
    
  ## GCA, Generalized Concentration Addition  
    mixture.result <- optimize(f = GCA.obj, interval = c(10^-20,1),
                               # Ci= rep(X[i],41),
                               Ci = Ci,
                               alpha = param[,1],
                               beta = param[,2],
                               model = antibiotox$sgl$model)


  GCA.val[i] <- mixture.result$minimum

  # End GCA
  
  ## Average prediction (assuming all are Weibull models)
  AVG.val[i] <- Weibull(log10(conc.total.i),alpha = mean(param[,1]),beta = mean(param[,2]))$ECmix

  ## IA, Independent Action
  
  IA.val[i] <- IA.pred(log10(Ci),alpha = param[,1],beta = param[,2])
}

df <- data.frame("X" = x.total, "CA" = CA.val,"GCA" = GCA.val,"AVG" = AVG.val,"IA" = IA.val)

df.single <- list()

for (i in 1:7){
  df.single[[i]] <- data.frame("X" = antibiotox[[i]]$x, "Y" = antibiotox[[i]]$y[,1])
}
df.single.melt <- melt(df.single,measure.vars = "Y")

df.ind <- df.single.melt[,c(1,3,4)]
colnames(df.ind) <- c("X","Y","variable")

df.melt <- melt(df,measure.vars = c("CA","GCA","AVG","IA"))
colnames(df.melt) <- c("X","variable","Y")
df.final <- rbind(df.ind,df.melt)

'%!in%' <- function(x,y)!('%in%'(x,y))

ggplot()+geom_point(data = df.final[df.final$variable %in% c("CA","GCA","AVG","IA"),],aes(X,Y,color = variable,shape = variable),size = 3)+
  geom_line(data = df.final[df.final$variable %in% c("CA","GCA","AVG","IA"),],aes(X,Y,color = variable),size = 1.5)+
  geom_point(data = df.final[df.final$variable %!in% c("CA","IA","GCA","AVG"),],aes(X,Y,group = variable))+
  geom_line(data = df.final[df.final$variable %!in% c("CA","IA","GCA","AVG"),],aes(X,Y,group = variable))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_color_viridis_d()+
  xlab("Concentration (mol/L)")+
  ylab("Effect")  



 
