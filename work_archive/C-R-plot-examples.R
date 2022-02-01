# Plotting the concentration-response curves of our data
library(ggplot2)
library(scales)
library(reshape2)

# The concentration-response uses log10 concentration
# We use the dataframe in the function to include the corresponding regular 
# space X
# We plot the x-axis using a log10 scaling 
tcplHillVal <- function(logc, tp, ga, gw, bt = 0) {
  
  bt + (tp - bt)/(1 + 10^((ga - logc)*gw))
  
}

## tcpl 3 parameter hill function inverse
tcplHillConc <- function(val, tp, ga, gw, bt = 0) {
  
  ga - log10((tp - bt)/(val - bt) - 1)/gw
  
}

# tp.mean <- get(load("/Volumes/SHAG/GeoTox/data/archived_outputs/tp_mean_example.RData"))
# AC50.mean <- get(load("/Volumes/SHAG/GeoTox/data/AC50_mean_example.RData"))
# slope.mean <- get(load("/Volumes/SHAG/GeoTox/data/slope_mean_example.RData"))

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
# If we want to include SD somehow in this example
# tp.sd <- get(load("/Volumes/SHAG/GeoTox/data/tp_sd_example.RData"))
# AC50.sd <- get(load("/Volumes/SHAG/GeoTox/data/AC50_sd_example.RData"))
# slope.sd <- get(load("/Volumes/SHAG/GeoTox/data/slope_sd_example.RData"))

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
