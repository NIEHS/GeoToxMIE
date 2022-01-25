

library(httk)
library(truncnorm)
library(tidyverse)

source("tcplHillVal.R", echo=FALSE)
# load data 
css.by.county <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/css_by_county_20220124.RData"))
cyp1a1_up.by.county <- get(load("/Volumes/SHAG/GeoTox/data/CYP1A1_by_county_20220124.RData"))
uchems <- get(load("/Volumes/SHAG/GeoTox/data/uchems_20220124.RData"))
inhalation.dose.by.county <- get(load("/Volumes/SHAG/GeoTox/data/inhalation_dose_by_county_20220124.RData"))

# chemical list
in.chems <-  c("98-86-2","92-87-5","92-52-4","117-81-7","133-06-2","532-27-4","133-90-4","57-74-9","510-15-6","94-75-7" ,
               "64-67-5","132-64-9","106-46-7","111-44-4","79-44-7","131-11-3","77-78-1","119-90-4","121-14-2","534-52-1", 
               "51-28-5","121-69-7","107-21-1","51-79-6","76-44-8","822-06-0","77-47-4","123-31-9","72-43-5" , 
               "101-77-9","56-38-2","82-68-8","87-86-5","1120-71-4", "114-26-1","91-22-5","96-09-3","95-80-7","584-84-9" ,
               "95-95-4","1582-09-8")

idx.chems <- uchems %in% in.chems
MC.iter <- 10^3
# drop the chemical column from the inhalation dose and cyp1a1 by county
for (i in 1:length(inhalation.dose.by.county)){
  inhalation.dose.by.county[[i]] <- inhalation.dose.by.county[[i]][,idx.chems]
  cyp1a1_up.by.county[[i]] <- cyp1a1_up.by.county[[i]][idx.chems,]
}

# Calculate the in-vitro dose using the Css
invitro.fun <- function(x){
  
  invitro <- inhalation.dose.by.county[[x]] * css.by.county[[x]]
  return(invitro)
}

invitro.conc.by.county <- lapply(1:length(inhalation.dose.by.county),invitro.fun)

# Proportion of each chemical by county

# proportion.fun <- function(x){
#   
#   proportion <- cyp1a1_up.by.county[[x]]$concentration_mean / sum(cyp1a1_up.by.county[[x]]$concentration_mean)
#   df  <- cbind(cyp1a1_up.by.county[[x]],"Propoprtion" = proportion)
#   return(df)
# }
# 
# cyp1a1_up.by.county <- lapply(1:length(cyp1a1_up.by.county),proportion.fun)

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
  tp.sd <- cyp1a1_up.by.county[[x]]$hill_tp_sd %>% as.numeric()
  AC50.mean <- cyp1a1_up.by.county[[x]]$hill_ga %>% as.numeric()
  AC50.sd <- cyp1a1_up.by.county[[x]]$hill_ga_sd %>% as.numeric()
  slope.mean <- cyp1a1_up.by.county[[x]]$hill_gw %>% as.numeric()
  slope.sd <- cyp1a1_up.by.county[[x]]$hill_gw_sd %>% as.numeric()
  
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


save(final.response.by.county,file = "/Volumes/SHAG/GeoTox/data/final_response_by_county_20220124.RData")


## Summary Stats

################################################################################################3
# calculate summary statistics from monte carlo
summary(unlist(final.response.by.county))

dr.median <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) median(x$DR))))
colnames(dr.median) <- "DR.median"
dr.mean <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) mean(x$DR))))
colnames(dr.mean) <- "DR.mean"
dr.95.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$DR, 0.95, na.rm = TRUE))))
colnames(dr.95.quantile) <- "DR.95.quantile"
dr.5.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$DR, 0.05,na.rm = TRUE))))
colnames(dr.5.quantile) <- "DR.5.quantile"

hq.median <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) median(x$HQ))))
colnames(hq.median) <- "HQ.median"
hq.mean <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) mean(x$HQ))))
colnames(hq.mean) <- "HQ.mean"
hq.95.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$HQ, 0.95, na.rm = TRUE))))
colnames(hq.95.quantile) <- "HQ.95.quantile"
hq.5.quantile <- as.data.frame(unlist(lapply(final.response.by.county, FUN = function(x) quantile(x$HQ, 0.05, na.rm = TRUE))))
colnames(hq.5.quantile) <- "HQ.5.quantile"

load ("/Volumes/SHAG/GeoTox/data/FIPS_by_county.RData")
FIPS <- as.data.frame(FIPS)

ivive.summary.df<- cbind(FIPS, dr.median, dr.mean, dr.95.quantile, dr.5.quantile,
                         hq.median, hq.mean, hq.95.quantile, hq.5.quantile)

write.csv(ivive.summary.df, "/Volumes/SHAG/GeoTox/data/mc_all_summary_df_20220124.csv")


# Plot the responses
convert.fun <- function(x){
  DR <- final.response.by.county[[x]]$DR
}

DR <- sapply(1:length(final.response.by.county),convert.fun)
DR.melt <- melt(DR)

ggplot(DR.melt[1:10^5,],aes(value,color = as.factor(Var2)))+geom_density()+
  theme(legend.position = "none")+scale_color_viridis_d()

convert.fun <- function(x){
  DR <- final.response.by.county[[x]]$HQ
}

HQ <- sapply(1:length(final.response.by.county),convert.fun)
HQ.melt <- melt(HQ)

ggplot(HQ.melt[1:100000,],aes(value,color = as.factor(Var2)))+geom_density()+
  theme(legend.position = "none")+scale_color_viridis_d()+scale_x_log10()


