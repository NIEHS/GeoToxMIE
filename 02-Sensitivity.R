##################################################################################################
# Sensitivity analysis - varying obesity 

# Libraries
library(tidyverse)
library(reshape2)
library(truncnorm)
library(sf)

# load data
load ("/Volumes/SHAG/GeoTox/data/FIPS_by_county.RData")
FIPS <- as.data.frame(FIPS)
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))


# helper code

local.path.functions <- c("/Volumes/messierkp/Projects/AEP-AOP/GeoToxMIE/helper_functions/")

# source("/Volumes/SHAG/GeoTox/R_functions/MC_pipeline/census.age.sim.R", echo=FALSE)
source(paste0(local.path.functions,"census-age-sim.R"), echo=FALSE)
source(paste0(local.path.functions,"sim-IR-BW.R"), echo=FALSE)
source(paste0(local.path.functions,"tcplHillVal.R"), echo=FALSE)
source(paste0(local.path.functions,"tcplHillVal_v2.R"), echo=FALSE)
source(paste0(local.path.functions,"tcplHillConc.R"), echo=FALSE)
source(paste0(local.path.functions,"tcplHillConc_v2.R"), echo=FALSE)
source(paste0(local.path.functions,"GCA-obj.R"),echo = FALSE)
source(paste0(local.path.functions,"IA-Pred.R"),echo = FALSE) 
source(paste0(local.path.functions,"ECmix-obj.R"),echo = FALSE) 

MC.iter <- 10^3
##################################################################################################
##### MC-ToxGeo-up-to-IVIVE ####

#Note: Run Order #2


# Load the dataframe with county FIPS, Pollutant Concentration, and EPA/ICE IVIVE data
county_cyp1a1_up <- get(load("/Volumes/SHAG/GeoTox/data/county_cyp1a1_up_20220201.RData"))

#### Convert the cyp1a1 data to a list by county ####
cyp1a1_up.by.county <- split(county_cyp1a1_up,as.factor(county_cyp1a1_up$FIPS))



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

# age.by.county.median<- lapply(age.by.county, FUN= function (x) median(x))
age.by.county.median <- lapply(age.by.county.median, FUN = function (x) replicate(MC.iter, x))


#### Inhalation Rate per body weight-by age #####
IR.by.county <- sim.IR.BW(MC.iter,age.by.county.median)


##### Simulate exposure concentrations ####
sim.chem.fun <- function(x){
  val <- matrix(0,nrow = MC.iter, ncol = nrow(cyp1a1_up.by.county[[x]]))
  # nrow(cyp1a1_up.by.county[[x]]) == 41 for this study
  print(x)
  if (nrow(cyp1a1_up.by.county[[x]])==0){
    
  }else{
    for (i in 1:ncol(val)){
      
      # val[,i] <- rnorm(MC.iter,cyp1a1_up.by.county[[x]]$concentration_mean[i],
      #                  cyp1a1_up.by.county[[x]]$concentration_sd[i])
      
      mean.i <- cyp1a1_up.by.county[[x]]$concentration_mean[i]
      
      
      # external concentration is held constant
      sd.i = 0
      
      if (mean.i==0){
        next
      }else if(mean.i > 0 & is.na(sd.i)){
        sim.i <- rep(mean.i,MC.iter)
      }else{
        sim.i <- rtruncnorm(MC.iter,a = 0, b= Inf,mean = mean.i ,
                            sd = sd.i)
      }
      
      
      
      val[,i] <- sim.i
      
    }
  }
  
  # val <- val  * replicate(42,IR.by.county[[x]])
  return(val)
}


external.dose.by.county <- lapply(1:length(cyp1a1_up.by.county),sim.chem.fun)

convert.fun <- function(x){
  print(x)
  (external.dose.by.county[[x]]/1000) * replicate(ncol(external.dose.by.county[[x]]),IR.by.county[[x]])
}

inhalation.dose.by.county <- lapply(1:length(external.dose.by.county),convert.fun)


####################################################################################

css.sensitivity.obesity <- get(load("/Volumes/SHAG/GeoTox/data/css_by_county_sensitivity_obesity_20220201.RData"))


####################################################################################

invitro.fun <- function(x){
  
  invitro <- inhalation.dose.by.county[[x]] * css.sensitivity.obesity[[x]]
  return(invitro)
}

invitro.conc.by.county <- lapply(1:length(inhalation.dose.by.county),invitro.fun)





run.dr.fun <- function(x){
  print(x)
  tp.mean <- cyp1a1_up.by.county[[x]]$tp
  tp.sd <- 0
  AC50.mean <- cyp1a1_up.by.county[[x]]$logAC50 
  AC50.sd <- 0
  
  # Simulation constraints based on TCPL pipeline
  # For the log10-AC50, we use the minimum of the observed response from the assays
  # or the minimum of the predicted exposure (in-vitro)
  resp.min <- cyp1a1_up.by.county[[x]]$resp_min
  resp.max <- cyp1a1_up.by.county[[x]]$resp_max
  logc.min <- cyp1a1_up.by.county[[x]]$logc_min
  logc.max <- cyp1a1_up.by.county[[x]]$logc_max
  nl <- length(tp.mean)
  

  mixture.response.interval <- c(-50,50)
  
  # Simulate the D-R parameters from a truncated normal distribution
  # We have MC.iter simulations for each chemical
  
  tp.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=0,b = resp.max*1.5,tp.mean,tp.sd)}) %>% t()
  AC50.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=logc.min - 2,b = logc.max + 0.5,AC50.mean,AC50.sd)}) %>% t()
  
  
  Ci <- invitro.conc.by.county[[x]]
  
  GCA.eff <- IA.eff <- GCA.HQ.10 <-IA.HQ.10 <- rep(NA,MC.iter)
  for (iter in 1:MC.iter){
    
    Cij <- Ci[iter,]
    tp.ij <- tp.sim[iter,]
    AC50.ij <- 10^AC50.sim[iter,]
    mixture.result <- optimize(f = GCA.obj, interval = mixture.response.interval,
                               Ci = Cij,
                               tp = tp.ij,
                               AC50 = AC50.ij)
    # print(mixture.result)
    GCA.eff[iter] <- exp(mixture.result$minimum)
    
    IA.eff[iter] <- IA.pred(Cij,tp.ij,AC50.ij)
    
    
    
    
    
    # Estimate the maximal response level of the mixture
    max.result <- optimize(f = GCA.obj, interval = mixture.response.interval,
                           Ci = Cij * 10^14,
                           tp = tp.ij,
                           AC50 = AC50.ij)
    
    max.response <- exp(max.result$minimum)
    
    
    
    
    E10 <- max.response * 0.1
    # Solve for EC10/AC10
    
    
    
    EC10.result <- optimize(f = ECx.obj, interval = c(-1000,1000),
                            E = E10,
                            Ci = Cij,
                            tp = tp.ij,
                            AC50 = AC50.ij)
    
    EC10.GCA <- EC10.result$minimum
    # print(EC10.GCA,digits = 3)
    # Calculate the AC10
    # This accounts for the uncertainty since it is based 
    # the randomly sampled top of the curve and logAC50 parameter
    E10.by.chem <- tp.ij * 0.1
    AC10.ij <- tcplHillConc_v2(E10.by.chem,tp.ij,AC50.ij,
                               rep(1,length(tp.ij)))
    
    
    sCij <- sum(Cij)
    if (EC10.GCA>0){
      GCA.HQ.10[iter] <- sCij / EC10.GCA
    }
    IA.HQ.10[iter]  <- sum(Cij / AC10.ij)
    
    
    
  }
  
  
  df <- data.frame("GCA.Eff"= GCA.eff,"IA.eff" = IA.eff, 
                   "GCA.HQ.10" = GCA.HQ.10,"IA.HQ.10" = IA.HQ.10)
  return(df)  
  
}

# This should be a list by county, with MC.iter elements in each list entry
final.response.by.county <- lapply(1:length(cyp1a1_up.by.county),run.dr.fun)


save(final.response.by.county,file = "/Volumes/SHAG/GeoTox/data/sensitivity_results_02_obesity.RData")

