######################################################
# Additive Dose Model
# By: Kristin Eccles
# Date: Oct 22nd, 2021
# Updated, Run, 12/1/2021
# Written in R Version 4.0.2
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
ice_data <- get(load("/Volumes/SHAG/GeoTox/data/210105_ICE_cHTS_invitrodbv33.Rdata"))
epa_data<-list.files(path = "/Volumes/SHAG/GeoTox/data/INVITRODB_V3_3_LEVEL5/",
                     pattern = "*LTEA_200730.csv", 
                     full.names = T) %>% 
  map_df(~read_csv(., col_types = cols(.default = "c"))) 
nata_df<- read.csv("/Volumes/SHAG/GeoTox/data/2014_NATA_CONCS.csv")
nata_chemicals <- read.csv("/Volumes/SHAG/GeoTox/data/NATA_pollutant_names_casrn.csv")
county_2014 <-st_read("/Volumes/SHAG/GeoTox/data/cb_2014_us_county_5m/cb_2014_us_county_5m.shp")


# TOX21 Data
# ice_data <- get(load("210105_ICE_cHTS_invitrodbv33.Rdata"))

# Usually a best practice not to rename variable like this
# Rename a new variable and remove the old one if need be
ice_data <- subset(ice_data, new_hitc == 1)

epa_data$m4id <-as.numeric(epa_data$m4id )

ice_epa_df<- left_join(ice_data, 
  epa_data[c("m4id", "hill_tp", "hill_tp_sd", "hill_ga", 
             "hill_ga_sd", "hill_gw", "hill_gw_sd","resp_max","logc_min","logc_max")], by= "m4id", keep=FALSE)

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
# subset tox21 data
ice_cyp1a1_up <- subset(ice_epa_df, aenm == "LTEA_HepaRG_CYP1A1_up" )

county_cyp1a1_up <- left_join(county_stack, ice_cyp1a1_up, by=c("casrn" = "casn"), keep= FALSE)
# county_cyp1a1_up<- county_cyp1a1_up[!is.na(county_cyp1a1_up$ACC), ]
# county_cyp1a1_up<- county_cyp1a1_up[!is.na(county_cyp1a1_up$concentration_mean), ]

# Save a working Rdata dataframe
save(county_cyp1a1_up,file = "/Volumes/SHAG/GeoTox/data/county_cyp1a1_up_2021201.RData")

