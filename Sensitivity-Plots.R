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

# Test using all the MC data - not the summary stats
#Load data
all <- read.csv("/Volumes/SHAG/GeoTox/data/httk_IVIVE/mc_all_summary_df.csv")
age <- read.csv("/Volumes/SHAG/GeoTox/data/httk_IVIVE/age_ivive_sensitivity_summary_df.csv")
obesity <- read.csv("/Volumes/SHAG/GeoTox/data/httk_IVIVE/obesity_ivive_sensitivity_summary_df.csv")
external <- read.csv("/Volumes/SHAG/GeoTox/data/httk_IVIVE/external_conc_sensitivity_summary_df.csv")
conc.resp <- read.csv("/Volumes/SHAG/GeoTox/data/httk_IVIVE/age_ivive_sensitivity_summary_df.csv")
tk <- read.csv("/Volumes/SHAG/GeoTox/data/httk_IVIVE/httk_ivive_sensitivity_summary_df.csv")

age$variable <- "Age"
obesity$variable <- "Obesity"
external$variable <- "External-Concentration"
conc.resp$variable <- "Concentration-Response"
all$variable <-"All"
tk$variable <-"TK parameters"

#combine df
sensitivity.df <- rbind(age, obesity, external, conc.resp, all, tk)
# sensitivity.df$variable = with(sensitivity.df, reorder("All", "Age", "Concentration-Response", "External-Concentration",
#                                             "Obesity", "TK parameters"))
# sensitivity.df <- mutate(variable = factor(variable, levels=c("All", "Age", "Concentration-Response", "External Concentration",
#                                           "Obesity", "Plasma Concentration")))
#   

################################################################################
# Exploratory data analysis
summary(all)
summary(age)
summary(obesity)
summary(external)
summary(tk)

################################################################################
#### Ridgeline plot ####

conc.resp.plot <-ggplot(sensitivity.df, aes(x = DR.median, y = as.factor(variable), fill = as.factor(variable))) +
  geom_density_ridges(alpha=0.5) + coord_cartesian(xlim = c(1e-10,1e-4))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")
conc.resp.plot

conc.resp.plot <-ggplot(sensitivity.df, aes(x = DR.median, y = as.factor(variable), fill = as.factor(variable))) +
  geom_boxplot(outlier.shape =  NA)+ scale_x_log10()+coord_cartesian(xlim = c(1e-15,1e-2))+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")
conc.resp.plot

# save_plot("/Volumes/SHAG/GeoTox/data/httk_IVIVE/CR_ridgeline_plot.tif", dr.plot, width = 20, height = 20, dpi = 200)

hq.plot <- ggplot(sensitivity.df, aes(x = (HQ.median), y = variable, fill = variable)) +
  geom_density_ridges(alpha=0.5) + scale_x_log10()+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Hazard Quotient")+
  ylab("Varying Parameter")
hq.plot

hq.plot <- ggplot(sensitivity.df, aes(x = (HQ.median), y = variable, fill = variable)) +
  geom_boxplot(outlier.shape = NA)+ scale_x_log10()+
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Hazard Quotient")+
  ylab("Varying Parameter")
hq.plot



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
  geom_density_ridges(alpha=0.5) +coord_cartesian(xlim = c(1e-20,10))+
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

