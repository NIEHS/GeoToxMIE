######################################################
# By: Kyle Messier
# Date: Oct 22nd, 2021
# Edits: Kristin Eccles
# Updated, Run, 01/24/2022
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

# load data
load("/Volumes/SHAG/GeoTox/data/final_response_by_county_20220124.RData")

load ("/Volumes/SHAG/GeoTox/data/FIPS_by_county.RData")
FIPS <- as.data.frame(FIPS)
  
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

dr.median <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) median(x$DR, na.rm = TRUE))))
colnames(dr.median) <- "DR.median"
dr.mean <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) mean(x$DR, na.rm = TRUE))))
colnames(dr.mean) <- "DR.mean"
dr.95.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$DR, 0.95, na.rm = TRUE))))
colnames(dr.95.quantile) <- "DR.95.quantile"
dr.5.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$DR, 0.05, na.rm = TRUE))))
colnames(dr.5.quantile) <- "DR.5.quantile"

hq.median <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) median(x$HQ, na.rm = TRUE))))
colnames(hq.median) <- "HQ.median"
hq.mean <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) mean(x$HQ, na.rm = TRUE))))
colnames(hq.mean) <- "HQ.mean"
hq.95.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$HQ, 0.95,na.rm = TRUE ))))
colnames(hq.95.quantile) <- "HQ.95.quantile"
hq.5.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$HQ, 0.05, na.rm = TRUE))))
colnames(hq.5.quantile) <- "HQ.5.quantile"

ivive.summary.df<- cbind(FIPS, dr.median, dr.mean, dr.95.quantile, dr.5.quantile,
                         hq.median, hq.mean, hq.95.quantile, hq.5.quantile)

write.csv(ivive.summary.df, "/Volumes/SHAG/GeoTox/data/mc_all_summary_df_20220124.csv")

#make data spatial
ivive_county_cyp1a1_up_sp<- left_join(county_2014, ivive.summary.df, by=c("countyid" = "FIPS"), keep=FALSE)
ivive_county_cyp1a1_up_sf <-st_as_sf(ivive_county_cyp1a1_up_sp)

################################################################################################
#### Plots ####
dr_cyp1a1_up_median <- ggplot(data = ivive_county_cyp1a1_up_sf, aes(fill=DR.median)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

dr_cyp1a1_up_95q <- ggplot(data = ivive_county_cyp1a1_up_sf, aes(fill=DR.95.quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

dr_cyp1a1_up_5q <- ggplot(data = ivive_county_cyp1a1_up_sf, aes(fill=DR.5.quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Log2 Fold Change", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

#### HAZARD QUOTIENT ####

hq_cyp1a1_up_median <- ggplot(data = ivive_county_cyp1a1_up_sf, aes(fill=HQ.median)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Hazard Quotient", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

hq_cyp1a1_up_95q <- ggplot(data = ivive_county_cyp1a1_up_sf, aes(fill=HQ.95.quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Hazard Quotient", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

hq_cyp1a1_up_5q <- ggplot(data = ivive_county_cyp1a1_up_sf, aes(fill=HQ.5.quantile)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  labs(fill="Sum")+
  scale_fill_distiller(name="Hazard Quotient", palette = "YlGnBu", direction = 1) +
  geom_sf(data = states, fill = NA, size=0.15)+
  theme(text = element_text(size = 14)) 

#Compile
dr.hq.figurev=ggarrange(dr_cyp1a1_up_5q, 
                        hq_cyp1a1_up_5q,
                       dr_cyp1a1_up_median, 
                       hq_cyp1a1_up_median,
                       dr_cyp1a1_up_95q,
                       hq_cyp1a1_up_95q,
                  labels = c("(A) 5th Percentile ", "(D) 5th Percentile",
                           "(B) Median", "(E) Median",
                             "(C) 95th Percentile", "(F) 95th Percentile"),
                  vjust = 1,
                  hjust = -0.5,
                  align = "hv",
                  ncol = 2, nrow = 3,
                  common.legend = FALSE,
                  legend = "right")
save_plot("/Volumes/SHAG/GeoTox/data/plots/dr_hq_figurev_20220126.tif", dr.hq.figurev, width = 40, height = 30, dpi = 200)







