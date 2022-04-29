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
# load data
sensitivity.ext.conc <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_ext_conc_20220415.RData"))
sensitivity.httk <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_httk_20220415.RData"))
sensitivity.obesity <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_obesity_20220415.RData"))
sensitivity.age <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_age_20220415.RData"))
baseline <- get(load("/Volumes/SHAG/GeoTox/data/final_response_by_county_202200415.RData"))
sensitivity.conc.resp <- get(load("/Volumes/SHAG/GeoTox/data/sensitivity_results_conc_resp_20220415.RData"))


#### GCA ####
sensitivity.GCA.Eff <- NULL
for (x in 1:length(sensitivity.ext.conc)){
  print(x)
  CR <- cbind(sensitivity.ext.conc[[x]]$GCA.Eff,
              sensitivity.httk[[x]]$GCA.Eff,
              sensitivity.obesity[[x]]$GCA.Eff,
              sensitivity.age[[x]]$GCA.Eff,
              sensitivity.conc.resp[[x]]$GCA.Eff,
              baseline[[x]]$GCA.Eff)
  
  sensitivity.GCA.Eff <- rbind(sensitivity.GCA.Eff,CR)
  

}

colnames(sensitivity.GCA.Eff) <- c("External Concentration","Toxicokinetic Parameters",
                              "Obesity","Age","Concentration-Response",
                              "Baseline")


CR.melt <- melt(sensitivity.GCA.Eff)

CR.melt$Var2 <- factor(CR.melt$Var2, levels = c( "External Concentration","Toxicokinetic Parameters",
                                                 "Concentration-Response","Age","Obesity", "Baseline"))

conc.resp.plot.GCA <-ggplot(CR.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  stat_density_ridges( #bandwidth = 0.5
                      geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      quantiles = 4, quantile_lines = FALSE
                      )+ 
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Predicted Log2 Fold Change \nmRNA Expression CYP1A1")+
  ylab("Varying Parameter")+
  theme_minimal()+
  coord_cartesian(clip = "off")+
  theme(text = element_text(size = 14), legend.position="none", axis.text=element_text(size=14),
        axis.title=element_text(size=14)) 
conc.resp.plot.GCA
#save_plot("/Volumes/SHAG/GeoTox/data/plots/Sensitivity_GCA_Eff__20220415.tif", conc.resp.plot.GCA,dpi = 200)


#### IA model ####
sensitivity.IA.eff <- NULL
for (x in 1:length(sensitivity.ext.conc)){
  print(x)
  CR <- cbind(sensitivity.ext.conc[[x]]$IA.eff,sensitivity.httk[[x]]$IA.eff,
              sensitivity.obesity[[x]]$IA.eff,sensitivity.age[[x]]$IA.eff,
              sensitivity.conc.resp[[x]]$IA.eff,baseline[[x]]$IA.eff)
  
  sensitivity.IA.eff <- rbind(sensitivity.IA.eff,CR)
  
  
}

colnames(sensitivity.IA.eff) <- c("External Concentration","Toxicokinetic Parameters",
                                  "Obesity","Age","Concentration-Response",
                                  "Baseline")

CR.IA.melt <- melt(sensitivity.IA.eff)

CR.IA.melt$Var2 <- factor(CR.IA.melt$Var2, levels = c( "External Concentration","Toxicokinetic Parameters",
                                             "Concentration-Response","Age","Obesity", "Baseline"))

conc.resp.plot.IA <-ggplot(CR.IA.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  stat_density_ridges( #bandwidth = 0.5
    geom = "density_ridges_gradient", calc_ecdf = TRUE, 
    quantiles = 4, quantile_lines = FALSE
    )  + 
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Predicted Log2 Fold Change \nmRNA Expression CYP1A1")+
  ylab("")+
  theme_minimal()+
  coord_cartesian(clip = "off")+
  theme(text = element_text(size = 14), legend.position="none", axis.text=element_text(size=14),
        axis.title=element_text(size=14), 
        axis.text.y=element_blank()) 
conc.resp.plot.IA
#save_plot("/Volumes/SHAG/GeoTox/data/plots/Sensitivity_IA_Eff__20220415.tif", conc.resp.plot.IA,  dpi = 200)

#### Risk Quotient ####
sensitivity.IA.HQ.10 <- NULL
for (x in 1:length(sensitivity.ext.conc)){
  print(x)
  CR <- cbind(sensitivity.ext.conc[[x]]$IA.HQ.10,sensitivity.httk[[x]]$IA.HQ.10,
              sensitivity.obesity[[x]]$IA.HQ.10,sensitivity.age[[x]]$IA.HQ.10,
              sensitivity.conc.resp[[x]]$IA.HQ.10,baseline[[x]]$IA.HQ.10)
  
  sensitivity.IA.HQ.10 <- rbind(sensitivity.IA.HQ.10,CR)
  
  
}

colnames(sensitivity.IA.HQ.10) <- c("External Concentration","Toxicokinetic Parameters",
                                    "Obesity","Age","Concentration-Response",
                                    "Baseline")


HQ.IA.melt <- melt(sensitivity.IA.HQ.10)

HQ.IA.melt$Var2 <- factor(HQ.IA.melt$Var2, levels = c( "External Concentration","Toxicokinetic Parameters",
                                                     "Concentration-Response","Age","Obesity", "Baseline"))

HQ.plot.IA <-ggplot(HQ.IA.melt, aes(x = value, y = as.factor(Var2), fill = as.factor(Var2))) +
  stat_density_ridges( #bandwidth = 0.5
    geom = "density_ridges_gradient", calc_ecdf = TRUE, 
    quantiles = 4, quantile_lines = FALSE
  )  + 
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Meidan CYP1A1 Risk Quotient")+
  ylab(" ")+
  theme_minimal()+
  coord_cartesian(clip = "off")+
  theme(text = element_text(size = 14), legend.position="none", axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        axis.text.y=element_blank()) 
HQ.plot.IA
#save_plot("/Volumes/SHAG/GeoTox/data/plots/Sensitivity_IA_HQ_20220415.tif", HQ.plot.IA)


#### Combine Plots ####
composite=ggarrange(conc.resp.plot.GCA , conc.resp.plot.IA , HQ.plot.IA, 
              labels = c( "A", "B", "C"),
              vjust = 1,
              align = "h",
              #hjust = -0.5,
              ncol = 3, nrow = 1,
              widths = c(1, 1, 1), 
              font.label = list(size = 20, color = "black", face = "bold"),
              common.legend = FALSE)

save_plot("sensitivity_composite_20220428.tif", composite, width = 35, height = 20, dpi = 300)


