######################################################
# By: Kyle Messier
# Date: Oct 22nd, 2021
# Edits: Kristin Eccles
# Updated, Run, 01/24/2022
# Written in R Version 4.0.2
######################################################
# Fifth script in the main CYCP1A1 analysis pipeline
# This script calculates the hazard quotient and  efficacy 
# under CA and IA assumptions by county

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
hill2.fit <- get(load("/Volumes/SHAG/GeoTox/data/Hill_2param_model_fit.RData"))


MC.iter <- 10^3

# Fix a few NA's from Css
for (i in 1:length(css.by.county)){
  for (j in 1:41){
    idx <- is.na(css.by.county[[i]][,j])
    if (sum(idx)>0){
      css.mean <- mean(css.by.county[[i]][!idx,j])
      css.by.county[[i]][idx,j] <- css.mean
    }
  }
}

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

  # mixture.response.interval <- c(min(resp.min) - 0.5 ,max(resp.max) * 4)
  mixture.response.interval <- c(-50,50)
  
  # Simulate the D-R parameters from a truncated normal distribution
  # We have MC.iter simulations for each chemical
  
  tp.sim <- sapply(1:MC.iter,function(i){rtruncnorm(nl,a=0,b = resp.max*1.2,tp.mean,tp.sd)}) %>% t()
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

save(final.response.by.county,file = "/Volumes/SHAG/GeoTox/data/final_response_by_county_20220201.RData")





