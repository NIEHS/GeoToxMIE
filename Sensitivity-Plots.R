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
all <- read.csv("mc_ivive_summary_df.csv")
age <- read.csv("/Volumes/SHAG/GeoTox/data/httk_IVIVE/age_ivive_sensitivity_summary_df.csv")
obesity <- read.csv("obesity_ivive_summary_df2.csv")
external <- read.csv("external_ivive_summary_df2.csv")
dr <- read.csv("dr_ivive_summary_df2.csv")
css <- read.csv("css_ivive_summary_df.csv")

age$variable <- "Age"
obesity$variable <- "Obesity"
external$variable <- "External Concentration"
dr$variable <- "Concentration Response"
all$variable <-"All"
css$variable <-"Plasma Concentration"

#combine df
sensitivity.df <- rbind(age, obesity, external, dr, all, css)
sensitivity.df$variable = with(sensitivity.df, reorder("All", "Age", "Concentration-Response", "External Concentration",
                                            "Obesity", "Plasma Concentration"))
sensitivity.df <- mutate(variable = factor(variable, levels=c("All", "Age", "Concentration-Response", "External Concentration",
                                          "Obesity", "Plasma Concentration"))) %>%
  

################################################################################
# Exploratory data analysis
summary(all)
summary(age)
summary(obesity)
summary(external)
summary(css)

################################################################################
#### Ridgeline plot ####

dr.plot <-ggplot(sensitivity.df, aes(x = DR.median, y = variable, fill = variable)) +
  geom_density_ridges(alpha=0.5) +
  theme_ridges() + 
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Median Log2 Fold Change mRNA Expression CYP1A1")+
  ylab("Varying Parameter")
dr.plot
save_plot("dr_ridgeline_plot.tif", dr.plot, width = 20, height = 20, dpi = 200)

hq.plot <- ggplot(sensitivity.df, aes(x = log10(HQ.median), y = variable, fill = variable)) +
  geom_density_ridges(alpha=0.5) +
  theme_ridges() + 
  scale_fill_viridis_d(option = "C")+
  theme(legend.position = "none")+
  xlab("Meidan Hazard Quotient")+
  ylab("Varying Parameter")
hq.plot
