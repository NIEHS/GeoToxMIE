##################################################################################################

#### Simulatation by county ####
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

#### Set up #### 

# Libraries
library(dplyr)
library(purrr)
library(readr)
library(reshape2)
library(stringr)
library(truncnorm)
library(sf)
library(ggplot2)
library(tigris)
library(maps)
library(sjPlot)
library(ggpubr)

# load data
load ("/Volumes/SHAG/GeoTox/data/FIPS_by_county.RData")
FIPS <- as.data.frame(FIPS)
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))


# helper code
source("/Volumes/SHAG/GeoTox/R_functions/MC_pipeline/census.age.sim.R", echo=FALSE)
source("/Volumes/SHAG/GeoTox/R_functions/MC_pipeline/sim.IR.BW.R", echo=FALSE)
source("/Volumes/SHAG/GeoTox/R_functions/MC_pipeline/tcplHillVal.R", echo=FALSE)

MC.iter <- 10^3
##################################################################################################
##### MC-ToxGeo-up-to-IVIVE ####

#Note: Run Order #2

#### Section A - Load data and source local functions ####
###
##
#
# Load the dataframe with county FIPS, Pollutant Concentration, and EPA/ICE IVIVE data
county_cyp1a1_up <- get(load("/Volumes/SHAG/GeoTox/data/county_cyp1a1_up_20211109.RData"))


#### MC Age iterations by county ####
# Age stratification

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

# sensitivity - hold age contant
age.by.county.median<- lapply(age.by.county, FUN= function (x) round(median(x)))
age.by.county.median <- lapply(age.by.county.median, FUN = function (x) replicate(MC.iter, x))

#### Inhalation Rate per body weight-by age #####
IR.by.county <- sim.IR.BW(MC.iter,age.by.county.median)


#### Convert the cyp1a1 data to a list by county ####
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

nchems <- nrow(cyp1a1_up.by.county[[1]])
convert.fun <- function(x){
  print(x)
  (external.dose.by.county[[x]]/1000) * replicate(42,IR.by.county[[x]])
}

inhalation.dose.by.county <- lapply(1:length(external.dose.by.county),convert.fun)

uchems <- cyp1a1_up.by.county[[1]]$casrn %>% unique()


css.by.county <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/css_by_county_20211209.RData"))

# To hold CSS constant
css.by.county.median<- lapply(css.by.county, FUN= function (x) median(x))
css.by.county.median <- lapply(css.by.county.median, FUN = function (x) replicate(MC.iter, x))

####################################################################################
#### MC-ToxGeo-Run-Risk-Measure ####
# Notes: Run order #5

# chemical list
in.chems <-  c("98-86-2","92-87-5","92-52-4","117-81-7","133-06-2","532-27-4","133-90-4","57-74-9","510-15-6","94-75-7" ,
               "64-67-5","132-64-9","106-46-7","111-44-4","79-44-7","131-11-3","77-78-1","119-90-4","121-14-2","534-52-1", 
               "51-28-5","121-69-7","107-21-1","51-79-6","76-44-8","822-06-0","77-47-4","123-31-9","72-43-5" , 
               "101-77-9","56-38-2","82-68-8","87-86-5","1120-71-4", "114-26-1","91-22-5","96-09-3","95-80-7","584-84-9" ,
               "95-95-4","1582-09-8")

idx.chems <- uchems %in% in.chems

# drop the chemical column from the inhalation dose and cyp1a1 by county
for (i in 1:length(inhalation.dose.by.county)){
  inhalation.dose.by.county[[i]] <- inhalation.dose.by.county[[i]][,idx.chems]
  cyp1a1_up.by.county[[i]] <- cyp1a1_up.by.county[[i]][idx.chems,]
}

# Calculate the in-vitro dose using the Css
invitro.fun <- function(x){
  
  invitro <- inhalation.dose.by.county[[x]] * css.by.county.median[[x]]
  return(invitro)
}

invitro.conc.by.county <- lapply(1:length(inhalation.dose.by.county),invitro.fun)

# Proportion of each chemical by county
proportion.fun <- function(x){
  # Get each MC iterations invtro dose sum across all chemicals
  mc.sums <- rowSums(invitro.conc.by.county[[x]])
  # Divide the individual concentrations by the sum of chemicals- for each MC
  proportion <- sweep(invitro.conc.by.county[[x]],1,mc.sums,"/")
  return(proportion)
}

proportion.by.county <- lapply(1:length(invitro.conc.by.county),proportion.fun)

run.dr.fun <- function(x){
  print(x)
  tp.mean <- cyp1a1_up.by.county[[x]]$hill_tp %>% as.numeric()
  AC50.mean <- cyp1a1_up.by.county[[x]]$hill_ga %>% as.numeric()
  slope.mean <- cyp1a1_up.by.county[[x]]$hill_gw %>% as.numeric()
  
  #To Hold Constant for Sensitivity Analysis 
  slope.sd <- 0
  tp.sd <- 0
  AC50.sd <-0
  
  # AC50.sd <- cyp1a1_up.by.county[[x]]$hill_ga_sd %>% as.numeric()
  # slope.sd <- cyp1a1_up.by.county[[x]]$hill_gw_sd %>% as.numeric()
  # tp.sd <- cyp1a1_up.by.county[[x]]$hill_tp_sd %>% as.numeric()
  # 

  
  # Simulation constraints based on TCPL pipeline
  # For the log10-AC50, we use the minimum of the observed response from the assays
  # or the minimum of the predicted exposure (in-vitro)
  resp.max <- cyp1a1_up.by.county[[x]]$resp_max %>% as.numeric()
  logc.min <- cyp1a1_up.by.county[[x]]$logc_min %>% as.numeric() - 2
  min.exposure <- inhalation.dose.by.county[[x]] %>% min(na.rm = TRUE) - 2
  min.constraint <- min(c(logc.min,min.exposure))
  nl <- length(tp.mean)
  # fix NAs  
  tp.sd[is.na(tp.sd)] <- mean(tp.sd[!is.na(tp.sd)])
  AC50.sd[is.na(AC50.sd)] <- mean(AC50.sd[!is.na(AC50.sd)])
  slope.sd[is.na(slope.sd)] <- mean(slope.sd[!is.na(slope.sd)])
  
  # Simulate the D-R parameters from a truncated normal distribution
  # We have MC.iter simulations for each chemical
  
  tp.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=0,b = resp.max*1.2,tp.mean,tp.sd)}) %>% t()
  AC50.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=min.constraint,b = 2,AC50.mean,AC50.sd)}) %>% t()
  slope.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=0.3,b = 8,slope.mean,slope.sd)}) %>% t()
  
  # Add the doses by chemical for each MC.iter
  dose.sum <- log10(rowSums(invitro.conc.by.county[[x]]))
  # Do the additivity based on the  concentration weighting
  tp.val <- rowSums(tp.sim * proportion.by.county[[x]])
  AC50.val <- rowSums(AC50.sim * proportion.by.county[[x]])
  slope.val <- rowSums(slope.sim * proportion.by.county[[x]])
  
  # CALCULATE THE DOSE RESPONSE!
  dose.response <- tcplHillVal(dose.sum,tp.val,AC50.val,slope.val)
  
  hazard.quotient <- rowSums(invitro.conc.by.county[[x]] / 10^AC50.sim)
  
  

  # print(cbind(hazard.quotient,hazard.quotient.weighted))
  df <- data.frame("DR"= dose.response,"HQ" = hazard.quotient,
                   "TP" = tp.val,"AC50" = AC50.val,"Slope" = slope.val)

  return(df)
}

# This should be a list by county, with MC.iter elements in each list entry
final.response.by.county <- lapply(1:length(cyp1a1_up.by.county),run.dr.fun)


save(final.response.by.county,file = "/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_external_conc_explore_input.RData")

