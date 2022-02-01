# Fit concentration-response data to the 2-parameter hill function
# for all 41 chemicals


### Load in the in-vitro concentration response data from ICE
ice_data <- get(load("/Volumes/SHAG/GeoTox/data/LTEA_HepaRG_CYP1A1_up 41 chems for Kyle 220131.RData"))

# split the data by chemical 
ice.data.by.chem <- split(ice_data,as.factor(ice_data$casn))

chem.names <- sapply(1:length(ice.data.by.chem),function(x) {unique(ice.data.by.chem[[x]]$chnm)})
chem.casn <- sapply(1:length(ice.data.by.chem),function(x) {unique(ice.data.by.chem[[x]]$casn)})

source("/Volumes/SHAG/GeoTox/R_functions/AEP-AOP/tcpl_my_2hill_fit.R")
source("/Volumes/messierkp/R_functions/AEP-AOP/tcpl_2hill_obj.R")


# Preallocate the results dataframe
df.results <- data.frame("tp" = rep(NA,length(ice.data.by.chem)),
                         "tp.sd" = rep(NA,length(ice.data.by.chem)),
                         "logAC50" = rep(NA,length(ice.data.by.chem)),
                         "logAC50.sd" = rep(NA,length(ice.data.by.chem)))
rownames(df.results) <- chem.names

for (i in 1:length(ice.data.by.chem)){
  print(i)
  my.2hill.mdl <- my_2hill_fit(ice.data.by.chem[[i]]$logc,ice.data.by.chem[[i]]$resp,log = TRUE)
  df.results$tp[i] <- my.2hill.mdl$par[1]
  df.results$tp.sd[i] <- my.2hill.mdl$sds[1]
  df.results$logAC50[i] <- my.2hill.mdl$par[2]
  df.results$logAC50.sd[i] <- my.2hill.mdl$sds[2]
}

df.results$casn <- chem.casn

# 5 tp and 1 log-AC50 sd are NA, due to the constraints and being near them
# We change those to equal the mean (like a Poisson mean )
df.results$tp.sd[is.na(df.results$tp.sd)] <- df.results$tp[is.na(df.results$tp.sd)]
df.results$logAC50.sd[is.na(df.results$logAC50.sd)] <- df.results$logAC50[is.na(df.results$logAC50.sd)]

save(df.results, file = "/Volumes/SHAG/GeoTox/data/Hill_2param_model_fit.RData")
#### Make some  plots to check

X <- 10^seq(log10(10^-3),log10(10^3),length.out = 10^3)
logX <- log10(X)

val.2hill <- matrix(NA,nrow = length(X),ncol = 41)
for (i in 1:41){
  val.2hill[,i] <- df.results$tp[i] /(1 + 10^( (df.results$logAC50[i] - logX) ))
}

colnames(val.2hill) <- chem.casn

val.comb <- data.frame("X" = X,val.2hill)


df.plot <-  reshape2::melt(val.comb,id.vars = "X")

ggplot(df.plot,aes(X,value,color = variable)) + geom_line()+ 
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
