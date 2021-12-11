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
sensitivity.df$variable = with(sensitivity.df, reorder("All", "Age", "Concentration-Response", "External-Concentration",
                                            "Obesity", "TK parameters"))
sensitivity.df <- mutate(variable = factor(variable, levels=c("All", "Age", "Concentration-Response", "External Concentration",
                                          "Obesity", "Plasma Concentration")))
  

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

hq.plot <- ggplot(sensitivity.df, aes(x = log10(HQ.median), y = variable, fill = variable)) +
  geom_density_ridges(alpha=0.5) +
  theme_ridges() + 
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Meidan Hazard Quotient")+
  ylab("Varying Parameter")
hq.plot
