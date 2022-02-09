#load libraries
library(psych)
library(ggplot2)
library(ggridges)
library(magrittr)
library(reshape2)
library(viridis)
library(plyr)
library(purrr)
library(dplyr)
library(scales)


## Sensitivity plots with each of the MC.iterations

external.concentration <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_results_external_conc.RData"))
sensitivity.httk <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_results_httk.RData"))
sensitivity.obesity <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_results_obesity.RData"))
sensitivity.age <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_results_age.RData"))
baseline <- get(load("/Volumes/SHAG/GeoTox/data/final_response_by_county_20211209.RData"))
sensitivity.conc.resp <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_results_conc_response.RData"))

sensitivity.CR <- NULL
for (x in 1:length(external.concentration)){
  print(x)
  CR <- cbind(external.concentration[[x]]$DR,sensitivity.httk[[x]]$DR,
              sensitivity.obesity[[x]]$DR,sensitivity.age[[x]]$DR,
              sensitivity.conc.resp[[x]]$DR,baseline[[x]]$DR)
  
  sensitivity.CR <- rbind(sensitivity.CR,CR)
  

}

colnames(sensitivity.CR) <- c("External-Concentration","TK-params",
                              "Obesity","Age","Concentration-Response",
                              "Baseline")


CR.melt <- melt(sensitivity.CR)

conc.resp.plot <-ggplot(CR.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  stat_density_ridges(bandwidth = 0.5,
                      geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = 4, quantile_lines = TRUE
  )  + 
  coord_cartesian(xlim = c(1e-20,10))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")
conc.resp.plot

conc.resp.plot <-ggplot(CR.melt, aes(x = value, fill = as.factor(Var2))) +
  geom_boxplot() +coord_cartesian(xlim = c(1e-20,10))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")
conc.resp.plot


conc.resp.facet <-ggplot(CR.melt, aes(value)) +
  geom_density()+scale_x_log10()+
  facet_wrap(~as.factor(Var2),scales = "free")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")
conc.resp.facet


sensitivity.HQ <- NULL
for (x in 1:length(external.concentration)){
  
  print(x)
  HQ <- cbind(external.concentration[[x]]$HQw,sensitivity.httk[[x]]$HQ,
              sensitivity.obesity[[x]]$HQ,sensitivity.age[[x]]$HQ,
              sensitivity.conc.resp[[x]]$HQ,baseline[[x]]$HQ)

  sensitivity.HQ <- rbind(sensitivity.HQ,HQ)
  
}

colnames(sensitivity.HQ) <- c("External-Concentration","TK-params",
                              "Obesity","Age","Concentration-Response",
                              "Baseline")

HQ.melt <- melt(sensitivity.HQ)


HQ.plot <-ggplot(HQ.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  geom_density_ridges(alpha=0.5) + 
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  coord_cartesian(xlim = c(1e-9,100))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Hazard Quotient CYP1A1")+
  ylab("Varying Parameter")
HQ.plot




###

CR.summary.df <- apply(sensitivity.CR,2,quantile,c(0.025,0.1,0.25,0.5,0.75,0.90,0.975))

CR.summary.melt <- melt(CR.summary.df)

ggplot() + geom_crossbar(data = CR.summary.melt[CR.summary.melt$Var1 == "50%" | CR.summary.melt$Var1 == "25%" | CR.summary.melt$Var1 == "75%",],
                         aes(x = as.factor(Var2),y = value))



conc.resp.plot <-ggplot(CR.melt, aes(x = value,y = as.factor(Var2), fill = as.factor(Var2))) +
  geom_boxplot(outlier.shape = NA) +coord_cartesian(xlim = c(1e-20,10))+
   scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")
conc.resp.plot


conc.resp.plot <-ggplot(HQ.melt, aes(x = value,y = as.factor(Var2), fill = as.factor(Var2))) +
  geom_boxplot(outlier.shape = NA) +coord_cartesian(xlim = c(1e-9,100))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Hazard Quotient of CYP1A1")+
  ylab("Varying Parameter")
conc.resp.plot

