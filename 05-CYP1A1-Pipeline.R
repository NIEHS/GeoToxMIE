######################################################
# By: Kyle Messier
# Date: Oct 22nd, 2021
# Edits: Kristin Eccles
# Updated, Run, 01/24/2022
# Written in R Version 4.0.2
######################################################
# Fifth script in the main CYCP1A1 analysis pipeline
# This script calculates the hazard quotient and GCA efficacy by county

# load libraries
library(httk)
library(truncnorm)
library(tidyverse)

#load code 
# specify the local path of the github repo
local.path.functions <- c("/Volumes/messierkp/Projects/AEP-AOP/GeoToxMIE/helper_functions/")
source(paste0(local.path.functions,"tcplHillVal.R"), echo=FALSE)
source(paste0(local.path.functions,"tcplHillVal_v2.R"), echo=FALSE)
source(paste0(local.path.functions,"tcplHillConc.R"), echo=FALSE)
source(paste0(local.path.functions,"tcplHillConc_v2.R"), echo=FALSE)
source(paste0(local.path.functions,"GCA-obj.R"),echo = FALSE)
source(paste0(local.path.functions,"IA-Pred.R"),echo = FALSE) 
source(paste0(local.path.functions,"ECmix-obj.R"),echo = FALSE) 


# load data 
css.by.county <- get(load("/Volumes/SHAG/GeoTox/data/css_by_county_20220201.RData"))
cyp1a1_up.by.county <- get(load("/Volumes/SHAG/GeoTox/data/CYP1A1_by_county_20220201.RData"))
inhalation.dose.by.county <- get(load("/Volumes/SHAG/GeoTox/data/inhalation_dose_by_county_20220201.RData"))


MC.iter <- 10^3



# Calculate the in-vitro dose using the Css
invitro.fun <- function(x){
  
  invitro <- inhalation.dose.by.county[[x]] * css.by.county[[x]]
  return(invitro)
}

invitro.conc.by.county <- lapply(1:length(inhalation.dose.by.county),invitro.fun)


run.dr.fun <- function(x){
  print(x)
  tp.mean <- cyp1a1_up.by.county[[x]]$tp
  tp.sd <- cyp1a1_up.by.county[[x]]$tp.sd
  AC50.mean <- cyp1a1_up.by.county[[x]]$logAC50 
  AC50.sd <- cyp1a1_up.by.county[[x]]$logAC50.sd

  
  # Simulation constraints based on TCPL pipeline
  # For the log10-AC50, we use the minimum of the observed response from the assays
  # or the minimum of the predicted exposure (in-vitro)
  resp.min <- cyp1a1_up.by.county[[x]]$resp_min
  resp.max <- cyp1a1_up.by.county[[x]]$resp_max
  logc.min <- cyp1a1_up.by.county[[x]]$logc_min
  logc.max <- cyp1a1_up.by.county[[x]]$logc_max
  nl <- length(tp.mean)

  mixture.response.interval <- c(min(resp.min) - 0.5 ,max(resp.max) * 4)
  
  # Simulate the D-R parameters from a truncated normal distribution
  # We have MC.iter simulations for each chemical
  
  tp.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=0,b = resp.max*1.2,tp.mean,tp.sd)}) %>% t()
  AC50.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=logc.min - 2,b = logc.max + 0.5,AC50.mean,AC50.sd)}) %>% t()

  
  Ci <- invitro.conc.by.county[[x]]

  
  GCA.eff <- IA.eff <- GCA.HQ.10 <- GCA.HQ.50 <- IA.HQ.10 <- IA.HQ.50 <- rep(NA,MC.iter)
  for (iter in 1:MC.iter){

    Cij <- Ci[iter,]
    tp.ij <- tp.sim[iter,]
    AC50.ij <- 10^AC50.sim[iter,]
    mixture.result <- optimize(f = GCA.obj, interval = mixture.response.interval,
                               Ci = Cij,
                               tp = tp.ij,
                               AC50 = AC50.ij)

    GCA.eff[iter] <- mixture.result$minimum
    
    IA.eff[iter] <- IA.pred(Cij,tp.ij,AC50.ij)
    
    
    # Estimate the maximal response level of the mixture
    max.result <- optimize(f = GCA.obj, interval = mixture.response.interval,
                               Ci = rep(10^3,length(tp.ij)),
                               tp = tp.ij,
                               AC50 = AC50.ij)
    
    max.response <- max.result$minimum
    
    E50 <- max.response * 0.5
    E10 <- max.response * 0.1
 
    # Solve for AC50 and AC10
      EC50.result <- optimize(f = ECx.obj, interval = c(0,1000),
                              E = E50,
                              Ci = Cij,
                              tp = tp.ij,
                              AC50 = AC50.ij)
      
      EC50.GCA<- EC50.result$minimum


    
      EC10.result <- optimize(f = ECx.obj, interval = c(0,1000),
                            E = E10,
                            Ci = Cij,
                            tp = tp.ij,
                            AC50 = AC50.ij)
      
      EC10.GCA <- EC10.result$minimum
    
      sCij <- sum(Cij)
      GCA.HQ.50[iter] <- sCij / EC50.GCA
      GCA.HQ.10[iter] <- sCij / EC10.GCA
      IA.HQ.50[iter]  <- sum(Cij)
    # Independent HQ - sum each HQ
    # CA HQ - calculate the ECx value - compare to sum_Ci
    
  }
  
  
  

  
  
  # Add the doses by chemical for each MC.iter
  dose.sum <- log10(rowSums(invitro.conc.by.county[[x]]))
  # Do the additivity based on the  concentration weighting
  tp.val <- rowSums(tp.sim * proportion.by.county[[x]])
  AC50.val <- rowSums(AC50.sim * proportion.by.county[[x]])

  # CALCULATE THE DOSE RESPONSE!
  dose.response <- tcplHillVal(dose.sum,tp.val,AC50.val,slope.val)
  
  hazard.quotient <- rowSums(invitro.conc.by.county[[x]] / 10^AC50.sim)

  df <- data.frame("DR"= dose.response,"HQ" = hazard.quotient)
  return(df)
}

# This should be a list by county, with MC.iter elements in each list entry
final.response.by.county <- lapply(1:length(cyp1a1_up.by.county),run.dr.fun)

save(final.response.by.county,file = "/Volumes/SHAG/GeoTox/data/final_response_by_county_20220201.RData")





