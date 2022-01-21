######################################################
# GeoTox Monte Carlo Analysis for CYP1A1
# By: Kristin Eccles
# Edits: Kyle P Messier
# Date: Oct 22nd, 2021
# Updated, Run, 01/12/2022
# QC by KPM, 12/7/2021
# Written in R Version 4.0.2
# Description: 
######################################################

# Libraries
library(ggplot2)
library(viridis)
library(dplyr)
library(tigris)
library(rgdal)
library(sf)
library(purrr)
library(readr)
library(reshape2)
library(sjPlot)
library(stringr)
library(viridisLite)
library(spdep)
library(maps)
library(ggpubr)
library(RColorBrewer)
library(httk)

# Import data

# From server (or local below)
#ice_data <- get(load("/Volumes/SHAG/GeoTox/data/210105_ICE_cHTS_invitrodbv33.Rdata"))
# epa_data<-list.files(path = "/Volumes/SHAG/GeoTox/data/INVITRODB_V3_3_LEVEL5/",
#                      pattern = "*LTEA_200730.csv", 
#                      full.names = T) %>% 
#   map_df(~read_csv(., col_types = cols(.default = "c"))) 

#Updated TK parameters
kdat_ice <- get(load("/Volumes/SHAG/GeoTox/data/220113_kdat_ice.RData"))
kdat_tcpl <- get(load("/Volumes/SHAG/GeoTox/data/220113_kdat_tcpl_details.RData"))

kdat_join <- left_join(kdat_tcpl, kdat_ice[,c("m4id", "new_ice_hitc")],by="m4id", keep=FALSE)

# limit to CYP1A1
kdat_ice_cyp<- subset(kdat_join, new_ice_hitc == 1 | new_ice_hitc == 3 )
kdat_ice_cyp <- subset(kdat_ice_cyp, aenm == "LTEA_HepaRG_CYP1A1_up")

# remove duplicates
kdat_ice_cyp<-unique(kdat_ice_cyp)
# rename to match variable names below
ice_epa_df <- kdat_ice_cyp

#epa_data$m4id <-as.numeric(epa_data$m4id )

# ice_epa_df<- left_join(ice_data, 
#                        epa_data[c("m4id", "hill_tp", "hill_tp_sd", "hill_ga", 
#                                   "hill_ga_sd", "hill_gw", "hill_gw_sd","resp_max","logc_min","logc_max")], by= "m4id", keep=FALSE)
# 

#limit to 
nata_df<- read.csv("/Volumes/SHAG/GeoTox/data/2014_NATA_CONCS.csv")
nata_chemicals <- read.csv("/Volumes/SHAG/GeoTox/data/NATA_pollutant_names_casrn.csv")
county_2014 <-st_read("/Volumes/SHAG/GeoTox/data/cb_2014_us_county_5m/cb_2014_us_county_5m.shp")

# There were some repeats of reading in data?!
#simplify names

nata_chemicals$web_name <- str_to_title(nata_chemicals$web_name, locale = "en")
nata_chemicals$web_name <- str_replace(nata_chemicals$web_name , " \\s*\\([^\\)]+\\)", "")

# aggregate to county
nata_county_mean <- aggregate(nata_df[,7:ncol(nata_df)] , by=list(nata_df$STCOFIPS) , FUN=mean) 
nata_county_sd <- aggregate(nata_df[,7:ncol(nata_df)] , by=list(nata_df$STCOFIPS) , FUN=sd) 

#melt
county_stack_mean <- melt(nata_county_mean, id=1)
colnames(county_stack_mean) <- cbind("FIPS", "chemical", "concentration_mean")

county_stack_sd <- melt(nata_county_sd, id=1)
colnames(county_stack_sd) <- cbind("FIPS", "chemical", "concentration_sd")

county_stack <- left_join(county_stack_mean, county_stack_sd, by=c("FIPS","chemical"), keep= FALSE)

# add casrn
county_stack <- left_join(county_stack, nata_chemicals, by=c("chemical" = "smoke_name"), keep= FALSE)

# Create a STATE ID 
county_stack$STATE <- county_stack$FIPS %>% as.character()
county_stack$STATE[nchar(county_stack$STATE)==4] <- str_pad(substr(county_stack$STATE[nchar(county_stack$STATE)==4],1,1),
                                                            width = 2,side = "left",pad = "0")
county_stack$STATE[nchar(county_stack$STATE)==5] <- substr(county_stack$STATE[nchar(county_stack$STATE)==5],1,2)


#limit to continental USA
county_stack<- subset(county_stack,  STATE != "02" )
county_stack<- subset(county_stack,  STATE != "15" )
county_stack<- subset(county_stack,  STATE != "60" )
county_stack<- subset(county_stack,  STATE != "66" )
county_stack<- subset(county_stack,  STATE != "69" )
county_stack<- subset(county_stack,  STATE != "72" )
county_stack<- subset(county_stack,  STATE != "78" )


################################################################################
#### CYP1A1 ####
county_cyp1a1_up <- left_join(county_stack, ice_epa_df, by=c("casrn" = "casn"), keep= FALSE)
#remove rows where NATA does not overlap with TOX21
#county_cyp1a1_up <- county_cyp1a1_up[!is.na(county_cyp1a1_up$aenm), ]
#county_cyp1a1_up <- county_cyp1a1_up[!is.na(county_cyp1a1_up$aenm), ]

# Save a working Rdata dataframe
save(county_cyp1a1_up,file = "/Volumes/SHAG/GeoTox/data/county_cyp1a1_up_20220121.RData")

