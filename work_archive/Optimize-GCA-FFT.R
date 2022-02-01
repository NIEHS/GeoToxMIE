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

tcplHillConc.Alternate <- function(e,tp,ga,gw){
  -1 * abs(ga * (tp/e -1))^(-1/gw)
  
}
# We

inv.tanh <- function(E,alpha,beta){
  atanh(E/alpha)/beta
}

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


adjust.obj <- function(param.new,param.old,x,w = rep(1,length(x))){
  # TRY THE SAME THING BUT WITH INVERSE FUNCTIONS INSTEAD
  # param.new = c(tp,AC50)
  # param.old = c(tp,AC50,slope)
  val.old <- tcplHillVal2(x,param.old[1],param.old[2],param.old[3])
  
  val.new <- tcplHillVal2(x,param.new[1],param.new[2],1)
  # val.new <- linear.forward(x,param.new[1],param.new[2])
  # val.new <- param.new[1] + param.new[2]*x + param.new[3]*x^2 + param.new[4]*x^3
  # val.new <- param.new[1]*tanh(param.new[2]*x)
  sum((w*(val.old - val.new))^2)
  
}

# Concentration addition objective function to find mixture effect
CA.obj <- function(ECmix,E,Ci,tp,AC50,slope){
  
  pi <- Ci / sum(Ci)
  ECi <- 10^tcplHillConc(E,tp,AC50,slope)
  # print(sum(is.na(ECi)))
  ca.val <- sum((pi * ECmix)/ECi,na.rm = )
  val <- (ca.val - 1)^2

  # print(c(ECmix,val))
  return(val)
  
}

GCA.obj <- function(E,Ci,tp,AC50,slope){
  
  # ECi <- tcplHillConc2(E,tp,AC50,slope)
  X <- 10^seq(log10(10^-3),log10(10^2),length.out = 10^3)
  
  for (i in 1:length(tp)){
    val <- tcplHillVal2(X,tp[i],AC50[i],slope[i])
    ECi <- fft(val,inverse = TRUE)
    
  }
  
  # ECi <- fftw::DCT(y) 
  # EC2 <- fftw::FFT(y)
  # EC3 <- fftw::IDCT(y)
  # EC4 <- fftw::IFFT(y)
   
  gca.val <- as.numeric(sum(Ci/ECi,na.rm = F))
  # print(sum(is.na(ECi)))
  val <- (gca.val - 1)^2
  # print(c("GCAval",gca.val,val))
   print(val)
  return(val)
  
}

IA.pred <- function(Ci,tp,AC50,slope){
  
  
  E <- tcplHillVal(Ci,tp,AC50,slope)
  
  1 - prod( (1 - E) )
  
}



##### END help functions #####

# Example hill model parameters (from the EPA/ICE database)
tp.mean <- c(1.75777662136449,
             2.64257335228456,
             1.2500912999803,
             3.17219027994445,
             3.44393885523638,
             0.956027091030637,
             2.45952587010046,
             5.23044416992704,
             8.75723605397531,
             4.43557862957007,
             1.70651562574913,
             2.19696943114655,
             0.970314753293581,
             1.09593968209227,
             1.62591751017549,
             0.952445370664443,
             5.14574689240698,
             9.70427307077555,
             4.69835271016736,
             4.90694556250327,
             4.35195253150607,
             2.37581751546163,
             1.4219266732615,
             1.87852653877426,
             3.81153247827415,
             2.35888545831149,
             1.47119714056452,
             9.05149720471746,
             3.04889555723495,
             4.81040784786756,
             5.55467567159707,
             5.50095390753827,
             3.24078787993538,
             1.32076005375849,
             1.13091863005194,
             2.34110862655768,
             1.99380635274463,
             7.35884622974272,
             3.20916264400714,
             3.33178272957434,
             5.66425708546636)

AC50.mean <- c(1.64296595037943,
               1.69511239992777,
               1.69746329917556,
               1.50608754296749,
               1.55075166332809,
               1.70507493706784,
               1.78931247162635,
               1.56357212970045,
               1.73191018044751,
               1.66873607917092,
               1.67757151234742,
               1.50662449422997,
               1.58076243146805,
               1.8261032207948,
               1.5330724108301,
               1.58153596317556,
               1.37690513671101,
               1.38742330175663,
               1.43065188607642,
               1.59868166946005,
               1.77383573921161,
               1.71423117136907,
               0.983974558803139,
               1.80652787865437,
               1.21005687314881,
               1.94026699327791,
               1.78011668273905,
               1.69369036758621,
               1.33801728153993,
               1.60246343717998,
               1.39135845813558,
               1.93154903774416,
               0.955440693012235,
               1.88822285530154,
               1.66855841664285,
               1.92180026051274,
               1.89280868461677,
               1.27197689658405,
               1.56394566837469,
               1.9324665453855,
               1.73903686680574)

slope.mean <- c(7.99969966935632,
                1.31329701493892,
                7.9758614678071,
                2.52388089970194,
                1.66773021828506,
                7.99876025022253,
                2.43874709812623,
                1.2496107837679,
                2.19158509004575,
                1.51872156597757,
                0.420310495512835,
                7.99987186843957,
                7.99999308064974,
                3.87183450239043,
                7.99999933890543,
                7.99999761703686,
                1.47318653372205,
                0.924477144252795,
                1.73715492680686,
                7.99996987820193,
                2.97637414092904,
                2.69805528421217,
                0.300000000113715,
                2.44899460904763,
                1.49361305266828,
                7.97752158888019,
                1.42482463931545,
                1.83508787957739,
                3.810875895109,
                1.29920040512761,
                1.17782930185412,
                7.99999995523108,
                1.25202823872039,
                5.09727088108956,
                7.99445094348028,
                7.99957905331678,
                2.04338236981753,
                0.808971549456102,
                1.37005005636836,
                7.99999999270583,
                1.63286331918594)


invitro.concentrations <- c(4.0686903230738e-06,
                            0,	
                            1.1411713731966e-07,
                            3.93549759097248e-07,	
                            1.0439212279716e-07,	
                            0,	
                            0,	
                            0,	
                            0,	
                            1.31104409034959e-05,	
                            0,	
                            1.27223297809921e-08,	
                            8.19608269084454e-07,	
                            0,	
                            0,	
                            6.59327077866334e-09,	
                            8.3482127306126e-12,	
                            0,	
                            7.27755146108946e-10,	
                            0,	
                            3.17085101920139e-11,	
                            6.77505096775645e-09,	
                            0.000195739321237987,	
                            0	,
                            0,	
                            0	,
                            3.65946983234194e-12,	
                            6.0233981058699e-08,	
                            0,	
                            0,	
                            0,	
                            1.86290746564753e-10,	
                            2.92381946690829e-09,	
                            0,	
                            0,	
                            0,	
                            0	,
                            0,	
                            1.02266414096264e-08,	
                            0	,
                            9.42741387002772e-07)


### Optimization appraoch
# Goal: Solve for E 
# 1) Initialize a value of E
# 2) Given the concentrations of compounds i = 1,..., n
# 3) Use inverse function to calculate ECxi 
# 4) Calculate the LHS of CA equation ( inv(sum(pi/ECxi)) = 1)
# 5) iterate until converges


X <- 10^seq(log10(10^-3),log10(10^2),length.out = 10^2)

init.temp <- matrix(NA,nrow = 41,ncol = length(X))
for (i in 1:length(X)){
  init.temp[,i] <- tcplHillVal2(X[i],tp.mean,10^AC50.mean,slope.mean)
}

chem.max <- apply(init.temp,1,max)
# ca.max <- min(chem.max)
ca.max <- max(chem.max)

# 
# E <- seq(10^-4,ca.max - 1e-4,length.out = 100)
# E.val <- rep(NA,length(X))
# assay.interval <- c(E[1],E[100]+1e-2)
 assay.interval <- c(10^-20,ca.max)
# assay.interval <- c(10^-20,10^3)
GCA.val <- rep(NA,length(X))
for (i in 1:length(X)){
  print(c(i))
  # Ci= sort(rnorm(41,mean = logX[i]))
  Ci <- rep(X[i],41)
  
  ## GCA, Generalized Concentration Addition  
  # GCA.obj <- function(E,Ci,tp,AC50,slope){
  mixture.result <- optimize(f = GCA.obj, interval = assay.interval,
                             # Ci= rep(X[i],41),
                             Ci = (Ci),
                             tp = tp.mean,
                             AC50 = 10^AC50.mean,
                             slope = slope.mean)
  
  
  GCA.val[i] <- mixture.result$minimum
 
}


# df1 <- data.frame("E" = E, "ECmix" = ECmix)
# df2 <- data.frame("E" = E.pred,"C" = conc.total.i$ECmix)
# ggplot()+geom_line(data = df1,aes(ECmix,E))+
#   geom_point(data = df1,aes(ECmix,E))+
#   geom_point(data = df2,aes(C,E),color = "red")

df1 <- data.frame("X" = X, "E" = GCA.val)
ggplot(df1,aes(X,E))+geom_line()+
  scale_x_log10()+ggtitle("GCA - slope substitution")
