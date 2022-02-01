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


## Begin the Dose-response modeling
x <- 1
  tp.mean <- cyp1a1_up.by.county[[x]]$hill_tp %>% as.numeric()
  AC50.mean <- cyp1a1_up.by.county[[x]]$hill_ga %>% as.numeric()
  slope.mean <- cyp1a1_up.by.county[[x]]$hill_gw %>% as.numeric()
  
  #To Hold Constant for Sensitivity Analysis 
  #slope.sd <- 0
  #tp.sd <- 0
  #AC50.sd <-0
  
  AC50.sd <- cyp1a1_up.by.county[[x]]$hill_ga_sd %>% as.numeric()
  slope.sd <- cyp1a1_up.by.county[[x]]$hill_gw_sd %>% as.numeric()
  tp.sd <- cyp1a1_up.by.county[[x]]$hill_tp_sd %>% as.numeric()
  
  # fix NAs  
  tp.sd[is.na(tp.sd)] <- mean(tp.sd[!is.na(tp.sd)])
  AC50.sd[is.na(AC50.sd)] <- mean(AC50.sd[!is.na(AC50.sd)])
  slope.sd[is.na(slope.sd)] <- mean(slope.sd[!is.na(slope.sd)])

  
  # Simulation constraints based on TCPL pipeline
  # For the log10-AC50, we use the minimum of the observed response from the assays
  # or the minimum of the predicted exposure (in-vitro)
  resp.max <- cyp1a1_up.by.county[[x]]$resp_max %>% as.numeric()
  logc.min <- cyp1a1_up.by.county[[x]]$logc_min %>% as.numeric() - 2
  min.constraint <- logc.min
  nl <- length(tp.mean)

  
  # Simulate the D-R parameters from a truncated normal distribution
  # We have MC.iter simulations for each chemical
  MC.iter <- 10^3
  tp.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=0,b = resp.max*1.2,tp.mean,tp.sd)}) %>% t()
  AC50.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=min.constraint,b = 2,AC50.mean,AC50.sd)}) %>% t()
  slope.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=0.3,b = 8,slope.mean,slope.sd)}) %>% t()
  
  
  ## Calculate and plot the distributions of each chemical
  logX <- log10(seq(10^-16,10^5,length.out = 1000))
  tp.val <- rowMeans(tp.sim)
  AC50.val <- rowMeans(AC50.sim)
  slope.val <- rowMeans(slope.sim)
  
  dose.response.DA <- NULL
  for (i in 1:length(logX)){
    DA <- tcplHillVal(logX[i],tp.val,AC50.val,slope.val)
    dose.response.DA <- rbind(dose.response.DA,DA)
    
  }
  
  CR.melt <- melt(dose.response.DA)
  conc.resp.plot <-ggplot(CR.melt, aes(x = value,y = as.factor(Var1))) +
    stat_density_ridges(calc_ecdf = TRUE, 
                        quantiles = 4, quantile_lines = TRUE,
                        geom = "density_ridges_gradient") +
    coord_cartesian(xlim = c(1e-2,10^5))+
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
    theme(legend.position = "none")+
    xlab("Median Log2 Fold Change mRNA Expression CYP1A1")
  conc.resp.plot
  
  
  
  # Add the doses by chemical for each MC.iter
  dose.sum <- log10(rowSums(invitro.conc.by.county[[x]]))
  # Do the additivity based on the  concentration weighting
  tp.val <- rowSums(tp.sim * proportion.by.county[[x]])
  AC50.val <- rowSums(AC50.sim * proportion.by.county[[x]])
  slope.val <- rowSums(slope.sim * proportion.by.county[[x]])
  
  # CALCULATE THE DOSE RESPONSE!
  # dose-addition
  dose.response.DA <- tcplHillVal(dose.sum,tp.val,AC50.val,slope.val)
  # response addition
  



# This should be a list by county, with MC.iter elements in each list entry
final.response.by.county <- lapply(1:length(cyp1a1_up.by.county),run.dr.fun)


save(final.response.by.county,file = "/Volumes/SHAG/GeoTox/data/httk_IVIVE/sensitivity_results_conc_response_RA.RData")

