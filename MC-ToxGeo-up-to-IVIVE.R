######################################################
# By: Kristin Eccles
# Date: Oct 22nd, 2021
# Edits: Kyle P Messier
# QC, 12/8/21, KPM
# Updated, Run, 01/12/2022
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
library(truncnorm)



# Load the dataframe with county FIPS, Pollutant Concentration, and EPA/ICE IVIVE data
county_cyp1a1_up <- get(load("/Volumes/SHAG/GeoTox/data/county_cyp1a1_up_20220121.RData"))

MC.iter <- 10^3

source("/Volumes/SHAG/GeoTox/R_functions/MC_pipeline/census.age.sim.R", echo=FALSE)
source("/Volumes/SHAG/GeoTox/R_functions/MC_pipeline/sim.IR.BW.R", echo=FALSE)
source("/Volumes/SHAG/GeoTox/R_functions/MC_pipeline/tcplHillVal.R", echo=FALSE)
source("/Volumes/SHAG/GeoTox/R_functions/MC_pipeline/Css.function.R", echo=FALSE)

#### MC Age iterations by county ####
#Age stratification
age.data <- read.csv("/Volumes/SHAG/GeoTox/data/cc-est2019-alldata.csv")
age.2014 <- subset(age.data, YEAR == 7) # 7/1/2014 Census population

# Create the county fips
age.county.fips <- str_pad(age.2014$COUNTY,width = 3,side = "left",pad = "0")
age.2014$FIPS <- paste0(age.2014$STATE,age.county.fips) %>% as.numeric()

# Census county name and FIPS change 
# https://www.ddorn.net/data/FIPS_County_Code_Changes.pdf
age.2014$FIPS[age.2014$FIPS==46102] <- 46113

census.age <- age.2014[,c("FIPS","AGEGRP","TOT_POP")]

idx.FIPS <- census.age$FIPS %in% county_cyp1a1_up$FIPS

census.age.overlap <- census.age[idx.FIPS,]

age.by.county <- census.age.sim(MC.iter,census.age.overlap)


#### Inhalation Rate per body weight-by age #####
IR.by.county <- sim.IR.BW(MC.iter,age.by.county)


#### CDC PLACES data on Obesity and Kidney disease prevalence ####

# Read in the Places data from the SHAG shared drive
places <- read.csv("/Volumes/SHAG/GeoTox/data/PLACES__County_Data__GIS_Friendly_Format___2020_release.csv")

# Calculate the standard deviation from the confidence internals
kidney.ci <- str_extract_all(places$KIDNEY_Crude95CI,pattern = "\\d+\\.\\d+") %>%  sapply(as.numeric) %>% t()
kidney.sd <- (kidney.ci[,2] - kidney.ci[,1])/3.92

obesity.ci <- str_extract_all(places$OBESITY_Crude95CI,pattern = "\\d+\\.\\d+") %>%  sapply(as.numeric) %>% t()
obesity.sd <- (obesity.ci[,2] - obesity.ci[,1])/3.92

places$KIDNEY_SD <- kidney.sd
places$OBESITY_SD <- obesity.sd

# Subset the columns we need
places <- places[,c("CountyFIPS","OBESITY_CrudePrev","OBESITY_SD")]

# Again, replacing the county FIPS because of a 2010 change
places$CountyFIPS[places$CountyFIPS==46102] <- 46113

# index and remove places counties outside conterminous US
idx.places <- places$CountyFIPS %in% unique(county_cyp1a1_up$FIPS)

places <- places[idx.places,]

# Simulate the Obesity data by county
obesity.binom.p <- lapply(1:length(places$OBESITY_CrudePrev),
                         function(x)rnorm(MC.iter,places$OBESITY_CrudePrev,places$OBESITY_SD)) 


obesity.by.county <- lapply(1:length(obesity.binom.p),
                           function(x)rbinom(MC.iter,size = 1,
                                             p = obesity.binom.p[[x]]/100)) 
# Restrict these to be only >=18
obesity.by.county <- lapply(1:length(obesity.by.county),function(x){
  val <- rep("Blank",length(obesity.by.county[[x]]))
  idx <- obesity.by.county[[x]]==0
  val[idx] = "Normal"
  idx <- obesity.by.county[[x]]==1
  val[idx] = "Obese"
  val[age.by.county[[x]]<18] = "Normal"
  return(val)
}
)


#### Convert the cyp1a1 data to a list by county ####
#county_cyp1a1_up$concentration_sd[is.na(county_cyp1a1_up$concentration_sd)] <- 0
cyp1a1_up.by.county <- split(county_cyp1a1_up,as.factor(county_cyp1a1_up$FIPS))

# Remove the chemicals that are not in the CYP1A1 AOP
NA.fun <- function(x){
  idx <- !is.na(cyp1a1_up.by.county[[x]]$hill_gw)
  val <- cyp1a1_up.by.county[[x]][idx,]
  return(val)
}
cyp1a1_up.by.county <- lapply(1:length(cyp1a1_up.by.county),NA.fun)


##### Simulate exposure concentrations ####
sim.chem.fun <- function(x){
  val <- matrix(0,nrow = MC.iter, ncol = 42)
  print(x)
  if (nrow(cyp1a1_up.by.county[[x]])==0){
    
  }else{
    for (i in 1:ncol(val)){
      
      # val[,i] <- rnorm(MC.iter,cyp1a1_up.by.county[[x]]$concentration_mean[i],
      #                  cyp1a1_up.by.county[[x]]$concentration_sd[i])
      
      mean.i <- cyp1a1_up.by.county[[x]]$concentration_mean[i]
      sd.i <- cyp1a1_up.by.county[[x]]$concentration_sd[i]
      
      if (mean.i==0){
        next
      }else if(mean.i > 0 & is.na(sd.i)){
        sim.i <- rep(mean.i,MC.iter)
      }else{
        sim.i <- rtruncnorm(MC.iter,a = 0, b= Inf,mean = mean.i,
                            sd = sd.i)
      }
      
      # sim.i[sim.i == "NaN"] <- 0
      val[,i] <- sim.i

    }
  }

  # val <- val  * replicate(42,IR.by.county[[x]])
  return(val)
}

# external.exposure.by.county <- lapply(1:length(cyp1a1_up.by.county),sim.chem.fun)
external.dose.by.county <- lapply(1:length(cyp1a1_up.by.county),sim.chem.fun)

#### Simulated parameters by county ####
# external.exposure.by.county
# age.by.county
# IR.by.county
# kidney.by.county
# obesity.by.county
# dose.response.by.county
# IVIVE.quant.by.county
# cyp1a1_up.by.county


#### air concentration from ug/m3 to mg/m3 then ####
# inhalation dose - multiply the air concentration by the inhalation rate per body weight 

nchems <- nrow(cyp1a1_up.by.county[[1]])
convert.fun <- function(x){
  print(x)
    (external.dose.by.county[[x]]/1000) * replicate(42,IR.by.county[[x]])
}

inhalation.dose.by.county <- lapply(1:length(external.dose.by.county),convert.fun)

uchems <- cyp1a1_up.by.county[[1]]$casrn %>% unique()

save(uchems,file = "/Volumes/SHAG/GeoTox/data/uchems_20211201.RData")
save(inhalation.dose.by.county,file = "/Volumes/SHAG/GeoTox/data/inhalation_dose_by_county_20220121.RData")
save(age.by.county,file = "/Volumes/SHAG/GeoTox/data/age_by_county_20220121.RData")
save(obesity.by.county,file = "/Volumes/SHAG/GeoTox/data/obesity_by_county_20220121.RData")
save(kidney.by.county,file = "/Volumes/SHAG/GeoTox/data/kidney_by_county_20220121.RData")
save(cyp1a1_up.by.county,file = "/Volumes/SHAG/GeoTox/data/CYP1A1_by_county_20220121.RData")
