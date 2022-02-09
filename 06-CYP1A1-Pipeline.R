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

# load data
load("/Volumes/SHAG/GeoTox/data/final_response_by_county_20220201.RData")

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
            mean = mean(value, na.rm = TRUE),
            x95_quantile = quantile(value, 0.95, na.rm = TRUE),
            x5_quantile = quantile(value, 0.05, na.rm = TRUE))%>%
            as.data.frame()

# add FIPS
ivive.summary.df <- left_join(ivive.summary.df, FIPS, by="id", keep=FALSE)
write.csv(ivive.summary.df, "/Volumes/SHAG/GeoTox/data/mc_all_summary_df_20220208.csv")

# make data spatial
ivive_county_cyp1a1_up_sp<- left_join(county_2014, ivive.summary.df, by=c("countyid" = "FIPS"), keep=FALSE)
ivive_county_cyp1a1_up_sf <-st_as_sf(ivive_county_cyp1a1_up_sp)

################################################################################################
#### Plots ####
#GCA.Eff, IA.eff, GCA.HQ.10, IA.HQ.10

#### GCA.Eff ####
GCA.Eff.median.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.Eff"), 
                              aes(fill=median)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

GCA.Eff.95.quantile.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.Eff"),
                                   aes(fill=x95_quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

GCA.Eff.5.quantile.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.Eff"),
                                  aes(fill=x5_quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

#Compile
GCA_Eff_figureh=ggarrange(GCA.Eff.5.quantile.plot, GCA.Eff.median.plot, GCA.Eff.95.quantile.plot,
                        labels = c("(A) 5th Percentile ", "(B) Median", "(C) 95th Percentile"),
                        vjust = 1,
                        hjust = -0.5,
                        align = "hv",
                        ncol = 3, nrow = 1,
                        common.legend = FALSE,
                        legend = "right")
save_plot("/Volumes/SHAG/GeoTox/data/plots/GCA_Eff_figureh_20220208.tif", GCA_Eff_figureh, width = 40, height = 7, dpi = 200)

#### IA.eff ####
IA.eff.median.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.eff"), 
                              aes(fill=median)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

IA.eff.95.quantile.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.eff"),
                                   aes(fill=x95_quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

IA.eff.5.quantile.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.eff"),
                                  aes(fill=x5_quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

#Compile
IA_eff_figureh=ggarrange(IA.eff.5.quantile.plot, IA.eff.median.plot, IA.eff.95.quantile.plot,
                          labels = c("(A) 5th Percentile ", "(B) Median", "(C) 95th Percentile"),
                          vjust = 1,
                          hjust = -0.5,
                          align = "hv",
                          ncol = 3, nrow = 1,
                          common.legend = FALSE,
                          legend = "right")
save_plot("/Volumes/SHAG/GeoTox/data/plots/IA_eff_figureh_20220208.tif", IA_eff_figureh, width = 40, height = 7, dpi = 200)

#### GCA.HQ.10 ####
GCA.HQ.10.median.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10"), 
                             aes(fill=median)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

GCA.HQ.10.95.quantile.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10"),
                                  aes(fill=x95_quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

GCA.HQ.10.5.quantile.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "GCA.HQ.10"),
                                 aes(fill=x5_quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

#Compile
GCA_HQ10_figureh=ggarrange(GCA.HQ.10.5.quantile.plot, GCA.HQ.10.median.plot, GCA.HQ.10.95.quantile.plot,
                         labels = c("(A) 5th Percentile ", "(B) Median", "(C) 95th Percentile"),
                         vjust = 1,
                         hjust = -0.5,
                         align = "hv",
                         ncol = 3, nrow = 1,
                         common.legend = FALSE,
                         legend = "right")
save_plot("/Volumes/SHAG/GeoTox/data/plots/GCA_HQ10_figureh_20220208.tif", GCA_HQ10_figureh, width = 40, height = 7, dpi = 200)


#### IA.HQ.10 ####
IA.HQ.10.median.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.HQ.10"), 
                                aes(fill=median)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

IA.HQ.10.95.quantile.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.HQ.10"),
                                     aes(fill=x95_quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

IA.HQ.10.5.quantile.plot <- ggplot(data = subset(ivive_county_cyp1a1_up_sf, health_measure == "IA.HQ.10"),
                                    aes(fill=x5_quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

#Compile
IA_HQ10_figureh=ggarrange(IA.HQ.10.5.quantile.plot, IA.HQ.10.median.plot, IA.HQ.10.95.quantile.plot,
                           labels = c("(A) 5th Percentile ", "(B) Median", "(C) 95th Percentile"),
                           vjust = 1,
                           hjust = -0.5,
                           align = "hv",
                           ncol = 3, nrow = 1,
                           common.legend = FALSE,
                           legend = "right")
save_plot("/Volumes/SHAG/GeoTox/data/plots/IA_HQ10_figureh_20220208.tif", IA_HQ10_figureh, width = 40, height = 7, dpi = 200)



