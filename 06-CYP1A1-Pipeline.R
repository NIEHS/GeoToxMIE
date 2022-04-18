######################################################
# By: Kristin Eccles
# Date: Oct 22nd, 2021
# Edits: Kyle Messier
# Updated, Run, 02/08/2022
# Written in R Version 4.0.2
######################################################

# load libraries
library(dplyr)
library(sf)
library(ggplot2)
library(tigris)
library(maps)
library(sjPlot)
library(viridisLite)
library(ggpubr)
library(reshape)
library(scales)

# load data
load("/Volumes/SHAG/GeoTox/data/final_response_by_county_202200415.RData")

load ("/Volumes/SHAG/GeoTox/data/FIPS_by_county.RData")
FIPS <- as.data.frame(FIPS)
# add list ids
FIPS$id <- 1:nrow(FIPS)

# Spatial Data
# state
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

#county 
county_2014 <-st_read("/Volumes/SHAG/GeoTox/data/cb_2014_us_county_5m/cb_2014_us_county_5m.shp")
county_2014$GEOID <- as.numeric(county_2014$GEOID)

#limit to continental USA
county_2014<- subset(county_2014,  STATEFP != "02" )
county_2014<- subset(county_2014,  STATEFP != "15" )
county_2014<- subset(county_2014,  STATEFP != "60" )
county_2014<- subset(county_2014,  STATEFP != "66" )
county_2014<- subset(county_2014,  STATEFP != "69" )
county_2014<- subset(county_2014,  STATEFP != "72" )
county_2014<- subset(county_2014,  STATEFP != "78" )
county_2014$countyid <-as.numeric(paste0(county_2014$STATEFP, county_2014$COUNTYFP))

################################################################################################
# calculate summary statistics from monte carlo
# unlist 
response.by.county <- melt(final.response.by.county)
colnames(response.by.county)<- cbind("health_measure", "value", "id")

#summarize
ivive.summary.df <- response.by.county %>%
  group_by(id, health_measure)%>%
  summarise(median = median(value , na.rm = TRUE),
            #mean = mean(value, na.rm = TRUE),
            x95_quantile = quantile(value, 0.95, na.rm = TRUE),
            x5_quantile = quantile(value, 0.05, na.rm = TRUE))%>%
            as.data.frame()
ivive.summary.df



mean(subset(ivive.summary.df, health_measure == "GCA.HQ.10")$median >1)*100
mean(subset(ivive.summary.df, health_measure == "GCA.HQ.10")$x5_quantile >1)*100
mean(subset(ivive.summary.df, health_measure == "GCA.HQ.10")$x95_quantile >1)*100

# add FIPS
ivive.summary.df <- left_join(ivive.summary.df, FIPS, by="id", keep=FALSE)
#write.csv(ivive.summary.df, "/Volumes/SHAG/GeoTox/data/mc_all_summary_df_20220415.csv")

# make data spatial
ivive.summary.df_stack <- melt(ivive.summary.df[,2:ncol(ivive.summary.df)], id.vars=c("FIPS", "health_measure"))
# order variable for plotting
ivive.summary.df_stack$variable = factor(ivive.summary.df_stack$variable, levels = c("x5_quantile","median","x95_quantile"))

# summary of health measure 
summary_stats <- ivive.summary.df_stack %>%
  group_by(variable, health_measure)%>%
  summarize(median=median(value),
            min=min(value), 
            max=max(value))%>%
  as.data.frame()
summary_stats
write.csv(summary_stats, "summary_stats.csv")


dose.labs <- c("GCA Response", "IA Response", "GCA RQ", "IA RQ")
names(dose.labs) <- c("GCA.Eff", "IA.eff", "GCA.HQ.10", "IA.HQ.10")

supp.labs <- c("5th Percentile", "Median", "95th Percentile")
names(supp.labs) <- c("x5_quantile", "median", "x95_quantile")

histogram_HM <- ggplot(data=ivive.summary.df_stack, aes(x=log10(value)))+
  geom_histogram(bins=500)+
  facet_grid(health_measure ~ variable,
             labeller=labeller(health_measure = dose.labs, variable =supp.labs))+
  theme_minimal()+
  ylab("Count")+
  xlab("Log10 Risk Metric Value")
histogram_HM
save_plot("/Volumes/SHAG/GeoTox/data/plots/health_metric_histogram_20220415.tif", histogram_HM, width = 30, height = 30, dpi = 200)




################################################################################################
#### Multi Panel Plots ####
#GCA.Eff, IA.eff, GCA.HQ.10, IA.HQ.10
ivive_county_cyp1a1_up_sp<- left_join(county_2014, ivive.summary.df_stack, by=c("countyid" = "FIPS"), keep=FALSE)
ivive_county_cyp1a1_up_sf <-st_as_sf(ivive_county_cyp1a1_up_sp)

plot.labs <- as_labeller(c(`x5_quantile` = "5th Percentile", `median` = "Median", `x95_quantile` = "95th Percentile"))

GCA.Eff.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.Eff"), 
                       aes(fill=value))+
  geom_sf(lwd = 0)+
  facet_wrap(~variable, nrow=2, ncol=3, labeller = plot.labs)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_viridis_c(name = "Response (Log2 FC)",direction = -1,option = "A",trans = "sqrt",
                       limits = c(0, 
                                  max(subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.Eff")$value)),
                       breaks = c(seq(min(subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.Eff")$value),
                                      max(subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.Eff")$value),by = 0.25)),
                       label = function(x) sprintf("%.2f", x))+
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 12), legend.text=element_text(size = 8)) 
GCA.Eff.plot
save_plot("/Volumes/SHAG/GeoTox/data/plots/GCA_Eff_figureh_20220415.tif", GCA.Eff.plot, width = 40, height = 7, dpi = 200)


IA.Eff.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.eff"), 
                       aes(fill=value))+
  geom_sf(lwd = 0)+
  facet_wrap(~variable, nrow=2, ncol=3, labeller = plot.labs)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_viridis_c(name = "Response (Log2 FC)",direction = -1,option = "A",trans = "sqrt",
                       limits = c(NA, 
                                  max(subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.eff")$value)),
                       breaks = c(seq(min(subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.eff")$value),
                                      max(subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.eff")$value),by = 0.25)),
                       label = function(x) sprintf("%.2f", x))+
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 12), legend.text=element_text(size = 8)) 
IA.Eff.plot
save_plot("/Volumes/SHAG/GeoTox/data/plots/IA_Eff_figureh_20220415.tif", IA.Eff.plot, width = 40, height = 7, dpi = 200)

GCA.HQ.10.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10"), 
                      aes(fill=value))+
  geom_sf(lwd = 0)+
  facet_wrap(~variable, nrow=2, ncol=3, labeller = plot.labs)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_viridis_c(name = "Risk Quotient",direction = -1,option = "A",trans = "sqrt",
                       limits = c(0, 
                                  round(max(subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10")$value))), 
                       breaks = c(round(min(subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10")$value)),
                                  10^seq(round(min(subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10")$value)),
                                         round(max(subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10")$value)),by = 0.5)), 
                       oob = scales::squish, 
                       label = function(x) sprintf("%.1f", x))+
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 12), legend.text=element_text(size = 8)) 
GCA.HQ.10.plot
save_plot("/Volumes/SHAG/GeoTox/data/plots/GCA_HQ10_figureh_20220415.tif", GCA.HQ.10.plot, width = 40, height = 7, dpi = 200)

HQ.1.10.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10" & value >1), 
                         aes(fill=value))+
  geom_sf(data = states, fill = "light grey", size=0.15)+
  geom_sf(lwd = 0)+
  geom_sf(fill=NA, lwd =0.05)+
  geom_sf(data = states, fill = NA, size=0.15)+
  facet_wrap(~variable, nrow=2, ncol=3, labeller = plot.labs)+
  theme_bw()+
  scale_fill_viridis_c(name = "RQ >1",direction = -1,option = "A")+
  theme(text = element_text(size = 12), legend.text=element_text(size = 8)) 
HQ.1.10.plot
save_plot("/Volumes/SHAG/GeoTox/data/plots/HQ10_1_figureh_20220415.tif", HQ.1.10.plot, width = 40, height = 7, dpi = 200)


IA.HQ.10.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.HQ.10"), 
                         aes(fill=value))+
  geom_sf(lwd = 0)+
  facet_wrap(~variable, nrow=2, ncol=3, labeller = plot.labs)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_viridis_c(name = "Risk Quotient",direction = -1,option = "A",trans = "sqrt",
                       limits = c(0, 
                                  round(max(subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.HQ.10")$value))), 
                       breaks = c(round(min(subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.HQ.10")$value)),
                                  10^seq(round(min(subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.HQ.10")$value)),
                                         round(max(subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.HQ.10")$value)),by = 0.5)), 
                       oob = scales::squish, 
                       label = function(x) sprintf("%.1f", x))+
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 12), legend.text=element_text(size = 8)) 
IA.HQ.10.plot
save_plot("/Volumes/SHAG/GeoTox/data/plots/IA_HQ10_figureh_20220415.tif", IA.HQ.10.plot, width = 40, height = 7, dpi = 200)

# 2x3 plot 
GCA=ggarrange(GCA.Eff.plot , GCA.HQ.10.plot, 
                     labels = c( "A", "B"),
                     vjust = 1,
                    align = "v",
                     #hjust = -0.5,
                     ncol = 1, nrow = 2,
                     font.label = list(size = 20, color = "black", face = "bold"),
                     common.legend = FALSE)
save_plot("/Volumes/SHAG/GeoTox/data/plots/GCA_composite_20220415.tif", GCA, width = 40, height = 15, dpi = 300)

IA=ggarrange(IA.Eff.plot , IA.HQ.10.plot, 
              labels = c( "A", "B"),
              vjust = 1,
              align = "v",
              #hjust = -0.5,
              ncol = 1, nrow = 2,
              font.label = list(size = 20, color = "black", face = "bold"),
              common.legend = FALSE)
save_plot("/Volumes/SHAG/GeoTox/data/plots/IA_composite_20220415.tif", IA, width = 40, height = 15, dpi = 300)

all=ggarrange(GCA.Eff.plot , GCA.HQ.10.plot, IA.Eff.plot , IA.HQ.10.plot, 
             labels = c( "(A) Concentration Addition - Response", "(B) Concentration Addition - Hazard Quotient",
                         "(C) Independent Action - Response", "(D) Indepdent Action - Harzard Quotient"),
             align = "v",
             #hjust = -0.5,
             ncol = 1, nrow = 4,
             font.label = list(size = 20, color = "black", face = "bold"),
             common.legend = FALSE)
save_plot("/Volumes/SHAG/GeoTox/data/plots/All_composite_20220415.tif", all, width = 40, height = 40, dpi = 300)


