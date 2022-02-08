# Fit concentration-response data to the 2-parameter hill function
# for all 41 chemicals


### Load in the in-vitro concentration response data from ICE
ice_data <- get(load("/Volumes/SHAG/GeoTox/data/LTEA_HepaRG_CYP1A1_up 41 chems for Kyle 220131.RData"))

# split the data by chemical 
ice.data.by.chem <- split(ice_data,as.factor(ice_data$casn))

chem.names <- sapply(1:length(ice.data.by.chem),function(x) {unique(ice.data.by.chem[[x]]$chnm)})
chem.casn <- sapply(1:length(ice.data.by.chem),function(x) {unique(ice.data.by.chem[[x]]$casn)})

source("/Volumes/messierkp/Projects/AEP-AOP/GeoToxMIE/helper_functions/tcpl_my_2hill_fit.R")
source("/Volumes/messierkp/Projects/AEP-AOP/GeoToxMIE/helper_functions/tcpl_2hill_obj.R")


# Preallocate the results dataframe
df.params <- data.frame("tp" = rep(NA,length(ice.data.by.chem)),
                         "tp.sd" = rep(NA,length(ice.data.by.chem)),
                         "logAC50" = rep(NA,length(ice.data.by.chem)),
                         "logAC50.sd" = rep(NA,length(ice.data.by.chem)),
                        "logc_min" = rep(NA,length(ice.data.by.chem)),
                        "logc_max" = rep(NA,length(ice.data.by.chem)),
                        "resp_min" = rep(NA,length(ice.data.by.chem)),
                        "resp_max" = rep(NA,length(ice.data.by.chem)),
                        "AIC" = rep(NA,length(ice.data.by.chem)))
rownames(df.params) <- chem.names

for (i in 1:length(ice.data.by.chem)){
  print(i)
  my.2hill.mdl <- my_2hill_fit(ice.data.by.chem[[i]]$logc,ice.data.by.chem[[i]]$resp,log = TRUE)
  df.params$tp[i] <- my.2hill.mdl$par[1]
  df.params$tp.sd[i] <- my.2hill.mdl$sds[1]
  df.params$logAC50[i] <- my.2hill.mdl$par[2]
  df.params$logAC50.sd[i] <- my.2hill.mdl$sds[2]
  df.params$logc_min[i] <- my.2hill.mdl$logc_min
  df.params$logc_max[i] <- my.2hill.mdl$logc_max
  df.params$resp_min[i] <- my.2hill.mdl$resp_min
  df.params$resp_max[i] <- my.2hill.mdl$resp_max
  df.params$AIC[i] <- my.2hill.mdl$AIC
}

df.params$casn <- chem.casn

# 5 tp and 1 log-AC50 sd are NA, due to the constraints and being near them
# We change those to equal the mean (like a Poisson mean )
df.params$tp.sd[is.na(df.params$tp.sd)] <- df.params$tp[is.na(df.params$tp.sd)]
df.params$logAC50.sd[is.na(df.params$logAC50.sd)] <- df.params$logAC50[is.na(df.params$logAC50.sd)]

save(df.params, file = "/Volumes/SHAG/GeoTox/data/Hill_2param_model_fit.RData")
#### Make some  plots to check

X <- 10^seq(log10(10^-3),log10(10^3),length.out = 10^3)
logX <- log10(X)

val.2hill <- matrix(NA,nrow = length(X),ncol = 41)
for (i in 1:41){
  val.2hill[,i] <- df.params$tp[i] /(1 + 10^( (df.params$logAC50[i] - logX) ))
}

colnames(val.2hill) <- chem.casn

val.comb <- data.frame("X" = X,val.2hill)


df.plot <-  reshape2::melt(val.comb,id.vars = "X")

ggplot(df.plot,aes(X,value,color = variable)) + geom_line()+ 
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))


for (i in 1:length(ice.data.by.chem)){
  df <- data.frame("X" = X,"Y" = val.2hill[,i])
  p <- ggplot()+ geom_line(data = df,aes(X,Y))+ 
    geom_point(data = ice.data.by.chem[[i]],aes(10^logc,resp),size = 2,color ="red")+ 
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
    ggtitle(ice.data.by.chem[[i]]$chnm[1])
  print(p)
  
}

