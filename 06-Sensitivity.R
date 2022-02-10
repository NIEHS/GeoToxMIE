#load libraries
library(psych)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(magrittr)
library(reshape2)
library(viridis)
library(plyr)
library(purrr)
library(dplyr)
library(scales)


## Sensitivity plots with each of the MC.iterations

sensitivity.ext.conc <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_05_ext_conc.RData"))
sensitivity.httk <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_03_httk.RData"))
sensitivity.obesity <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_02_obesity.RData"))
sensitivity.age <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_01_age.RData"))
baseline <- get(load("/Volumes/SHAG/GeoTox/data/final_response_by_county_20220201.RData"))
sensitivity.conc.resp <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_04_conc_resp.RData"))


sensitivity.GCA.Eff <- NULL
for (x in 1:length(sensitivity.ext.conc)){
  print(x)
  CR <- cbind(sensitivity.ext.conc[[x]]$GCA.Eff,sensitivity.httk[[x]]$GCA.Eff,
              sensitivity.obesity[[x]]$GCA.Eff,sensitivity.age[[x]]$GCA.Eff,
              sensitivity.conc.resp[[x]]$GCA.Eff,baseline[[x]]$GCA.Eff)
  
  sensitivity.GCA.Eff <- rbind(sensitivity.GCA.Eff,CR)
  

}

colnames(sensitivity.GCA.Eff) <- c("External-Concentration","TK-params",
                              "Obesity","Age","Concentration-Response",
                              "Baseline")


CR.melt <- melt(sensitivity.GCA.Eff)

conc.resp.plot.GCA <-ggplot(CR.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  stat_density_ridges( #bandwidth = 0.5
                      geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = 4, quantile_lines = TRUE
  )  + 
   coord_cartesian(xlim = c(10^-10,10^1))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")+ggtitle("Sensitivity Results - Efficacy - GCA model")
conc.resp.plot.GCA
save_plot("/Volumes/SHAG/GeoTox/data/plots/Sensitivity_GCA_Eff_20220209.tif", conc.resp.plot.GCA,dpi = 200)



sensitivity.IA.eff <- NULL
for (x in 1:length(sensitivity.ext.conc)){
  print(x)
  CR <- cbind(sensitivity.ext.conc[[x]]$IA.eff,sensitivity.httk[[x]]$IA.eff,
              sensitivity.obesity[[x]]$IA.eff,sensitivity.age[[x]]$IA.eff,
              sensitivity.conc.resp[[x]]$IA.eff,baseline[[x]]$IA.eff)
  
  sensitivity.IA.eff <- rbind(sensitivity.IA.eff,CR)
  
  
}

colnames(sensitivity.IA.eff) <- c("External-Concentration","TK-params",
                                   "Obesity","Age","Concentration-Response",
                                   "Baseline")


CR.IA.melt <- melt(sensitivity.IA.eff)

conc.resp.plot.IA <-ggplot(CR.IA.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  stat_density_ridges( #bandwidth = 0.5
    geom = "density_ridges_gradient", calc_ecdf = TRUE, 
    quantiles = 4, quantile_lines = TRUE
  )  + 
  coord_cartesian(xlim = c(10^-10,10^1))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")+ggtitle("Sensitivity Results - Efficacy - IA model")
conc.resp.plot.IA
save_plot("/Volumes/SHAG/GeoTox/data/plots/Sensitivity_IA_Eff_20220209.tif", conc.resp.plot.IA,  dpi = 200)

########## Hazard Quotients

sensitivity.GCA.HQ.10<- NULL
for (x in 1:length(sensitivity.ext.conc)){
  print(x)
  CR <- cbind(sensitivity.ext.conc[[x]]$GCA.HQ.10,sensitivity.httk[[x]]$GCA.HQ.10,
              sensitivity.obesity[[x]]$GCA.HQ.10,sensitivity.age[[x]]$GCA.HQ.10,
              sensitivity.conc.resp[[x]]$GCA.HQ.10,baseline[[x]]$GCA.HQ.10)
  
  sensitivity.GCA.HQ.10<- rbind(sensitivity.GCA.HQ.10,CR)
  
  
}

colnames(sensitivity.GCA.HQ.10) <- c("External-Concentration","TK-params",
                                   "Obesity","Age","Concentration-Response",
                                   "Baseline")


HQ.GCA.melt <- melt(sensitivity.GCA.HQ.10)

HQ.plot.GCA <-ggplot(HQ.GCA.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  stat_density_ridges( #bandwidth = 0.5
    geom = "density_ridges_gradient", calc_ecdf = TRUE, 
    quantiles = 4, quantile_lines = TRUE
  )  + 
  coord_cartesian(xlim = c(10^-10,10^1))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Hazard Quotient, CYP1A1")+
  ylab("Varying Parameter")+ggtitle("Sensitivity Results - HQ - GCA model")
HQ.plot.GCA
save_plot("/Volumes/SHAG/GeoTox/data/plots/Sensitivity_GCA_HQ_20220209.tif", HQ.plot.GCA, dpi = 200)



sensitivity.IA.HQ.10 <- NULL
for (x in 1:length(sensitivity.ext.conc)){
  print(x)
  CR <- cbind(sensitivity.ext.conc[[x]]$IA.HQ.10,sensitivity.httk[[x]]$IA.HQ.10,
              sensitivity.obesity[[x]]$IA.HQ.10,sensitivity.age[[x]]$IA.HQ.10,
              sensitivity.conc.resp[[x]]$IA.HQ.10,baseline[[x]]$IA.HQ.10)
  
  sensitivity.IA.HQ.10 <- rbind(sensitivity.IA.HQ.10,CR)
  
  
}

colnames(sensitivity.IA.HQ.10) <- c("External-Concentration","TK-params",
                                  "Obesity","Age","Concentration-Response",
                                  "Baseline")


HQ.IA.melt <- melt(sensitivity.IA.HQ.10)

HQ.plot.IA <-ggplot(HQ.IA.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  stat_density_ridges( #bandwidth = 0.5
    geom = "density_ridges_gradient", calc_ecdf = TRUE, 
    quantiles = 4, quantile_lines = TRUE
  )  + 
  coord_cartesian(xlim = c(10^-10,10^1))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Hazard Quotient CYP1A1")+
  ylab("Varying Parameter")+ggtitle("Sensitivity Results - HQ - IA model")
HQ.plot.IA
save_plot("/Volumes/SHAG/GeoTox/data/plots/Sensitivity_IA_HQ_20220209.png", HQ.plot.IA)



