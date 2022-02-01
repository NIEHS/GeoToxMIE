##################################################################################################
# Sensitivity analysis - do all of them together - change paramaters by
# 1, 2, and 3 standard deviations of their given distributions
# Written: KPM 12/13/21

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

# The standard normal quantiles
std.norm.quants <- pnorm(c(-3,-2,-1,0,1,2,3))

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

# SENSITIVITY AGE - uncomment hold age constant
age.by.county.quant1<- lapply(age.by.county, FUN= function (x) quantile(x,std.norm.quants[1]))
age.by.county.quant2<- lapply(age.by.county, FUN= function (x) quantile(x,std.norm.quants[2]))
age.by.county.quant3<- lapply(age.by.county, FUN= function (x) quantile(x,std.norm.quants[3]))
age.by.county.quant4<- lapply(age.by.county, FUN= function (x) quantile(x,std.norm.quants[4]))
age.by.county.quant5<- lapply(age.by.county, FUN= function (x) quantile(x,std.norm.quants[5]))
age.by.county.quant6<- lapply(age.by.county, FUN= function (x) quantile(x,std.norm.quants[6]))
age.by.county.quant7<- lapply(age.by.county, FUN= function (x) quantile(x,std.norm.quants[7]))

#### Inhalation Rate per body weight-by age #####
IR.by.county.quant1 <- sim.IR.BW(MC.iter,age.by.county.quant1)
IR.by.county.quant2 <- sim.IR.BW(MC.iter,age.by.county.quant2)
IR.by.county.quant3 <- sim.IR.BW(MC.iter,age.by.county.quant3)
IR.by.county.quant4 <- sim.IR.BW(MC.iter,age.by.county.quant4)
IR.by.county.quant5 <- sim.IR.BW(MC.iter,age.by.county.quant5)
IR.by.county.quant6 <- sim.IR.BW(MC.iter,age.by.county.quant6)
IR.by.county.quant7 <- sim.IR.BW(MC.iter,age.by.county.quant7)

#### CDC PLACES data on Obesity prevalence ####

# Read in the Places data from the SHAG shared drive
places <- read.csv("/Volumes/SHAG/GeoTox/data/PLACES__County_Data__GIS_Friendly_Format___2020_release.csv")

# Calculate the standard deviation from the confidence internals

obesity.ci <- str_extract_all(places$OBESITY_Crude95CI,pattern = "\\d+\\.\\d+") %>%  sapply(as.numeric) %>% t()
obesity.sd <- (obesity.ci[,2] - obesity.ci[,1])/3.92

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

#SENSITIVITY OBESITY - uncomment to hold contant
# obesity.by.county.median<- lapply(obesity.by.county, FUN= function (x) round(median(x)))
# obesity.by.county <- lapply(obesity.by.county.median, FUN = function (x) replicate(MC.iter, x))

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


#SENSITIVITY EXTERNAL CONCENTRATION- UNCOMMENT TO HOLD CONTSTANT
      #sd.i <- cyp1a1_up.by.county[[x]]$concentration_sd[i]
      sd.i = 0
      
      if (mean.i==0){
        next
      }else if(mean.i > 0 & is.na(sd.i)){
        sim.i <- rep(mean.i,MC.iter)
      }else{
        sim.i <- rtruncnorm(MC.iter,a = 0, b= Inf,mean = mean.i ,
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

####################################################################################

# To hold CSS constant
# load the sensitivity IVIVE data - css.sensitivity.age
load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/css_by_county_sensitivity_obesity.RData")


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
# Css pre-calculated to only include variability from age 
invitro.fun <- function(x){
  
  invitro <- inhalation.dose.by.county[[x]] * css.sensitivity.obesity[[x]]
  return(invitro)
}

invitro.conc.by.county <- lapply(1:length(inhalation.dose.by.county),invitro.fun)

# Proportion of each chemical by css

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
  
  #Sensitivity Analysis- Uncomment to hold constant 
  #AC50.sd <- cyp1a1_up.by.county[[x]]$hill_ga_sd %>% as.numeric()
  #slope.sd <- cyp1a1_up.by.county[[x]]$hill_gw_sd %>% as.numeric()
  #tp.sd <- cyp1a1_up.by.county[[x]]$hill_tp_sd %>% as.numeric()
  
  
  slope.sd <- 0
  tp.sd <- 0
  AC50.sd <-0
  
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
  
  df <- data.frame("DR"= dose.response,"HQ" = hazard.quotient)
  return(df)
}

# This should be a list by county, with MC.iter elements in each list entry
final.response.by.county <- lapply(1:length(cyp1a1_up.by.county),run.dr.fun)


save(final.response.by.county,file = "/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_results_obesity.RData")

######################################################################################
#load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_results_obesity.RData")
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

# calculate summary statistics from monte carlo

dr.median <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) median(x$DR))))
colnames(dr.median) <- "DR.median"
dr.mean <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) mean(x$DR))))
colnames(dr.mean) <- "DR.mean"
dr.95.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$DR, 0.95))))
colnames(dr.95.quantile) <- "DR.95.quantile"
dr.5.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$DR, 0.05))))
colnames(dr.5.quantile) <- "DR.5.quantile"

hq.median <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) median(x$HQ))))
colnames(hq.median) <- "HQ.median"
hq.mean <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) mean(x$HQ))))
colnames(hq.mean) <- "HQ.mean"
hq.95.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$HQ, 0.95))))
colnames(hq.95.quantile) <- "HQ.95.quantile"
hq.5.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$HQ, 0.05))))
colnames(hq.5.quantile) <- "HQ.5.quantile"

ivive.summary.df<- cbind(FIPS, dr.median, dr.mean, dr.95.quantile, dr.5.quantile,
                         hq.median, hq.mean, hq.95.quantile, hq.5.quantile)
summary(ivive.summary.df)
write.csv(ivive.summary.df, "/Volumes/SHAG/GeoTox/data/httk_IVIVE/obesity_ivive_sensitivity_summary_df.csv")
####################################################################################################
#### DOSE RESPONSE ####

states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

ivive_county_cyp1a1_up_sp<- left_join(county_2014, ivive.summary.df, by=c("countyid" = "FIPS"), keep=FALSE)
ivive_county_cyp1a1_up_sf <-st_as_sf(ivive_county_cyp1a1_up_sp)

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

##### Compile Figures ####

#Compile
sensitivity.age.dr.hq.figure = ggarrange(dr_cyp1a1_up_5q, hq_cyp1a1_up_5q,
                       dr_cyp1a1_up_median,hq_cyp1a1_up_median,
                       dr_cyp1a1_up_95q, hq_cyp1a1_up_95q,
                       labels = c("(A) 5th Percentile ", "(D) 5th Percentile",
                                  "(B) Median", "(E) Median",
                                  "(C) 95th Percentile", "(F) 95th Percentile"),
                       vjust = 1,
                       hjust = -0.5,
                       align = "hv",
                       ncol = 2, nrow = 3,
                       common.legend = FALSE,
                       legend = "right")

save_plot("/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity.age.dr.hq.figure2.tif", sensitivity.age.dr.hq.figure, width = 40, height = 30, dpi = 200)



