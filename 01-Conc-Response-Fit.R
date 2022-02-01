# Fit concentration-response data to the 2-parameter hill function
# Compare to the 3-parameter ICE data
# Check the standard deviation estimates, etc. 



# ICE-TCPL results to compare against
county_cyp1a1_up <- get(load("/Volumes/SHAG/GeoTox/data/archived_outputs/county_cyp1a1_up_20211109.RData"))

in_chems <-  c("98-86-2","92-87-5","92-52-4","117-81-7","133-06-2","532-27-4",
               "133-90-4","57-74-9","510-15-6","94-75-7","64-67-5","132-64-9",
               "106-46-7","111-44-4","79-44-7","131-11-3","77-78-1","119-90-4",
               "121-14-2","534-52-1","51-28-5","121-69-7","107-21-1","51-79-6",
               "76-44-8","822-06-0","77-47-4","123-31-9","72-43-5","101-77-9",
               "56-38-2","82-68-8","87-86-5","1120-71-4","114-26-1","91-22-5",
               "96-09-3","95-80-7","584-84-9","95-95-4","1582-09-8")

ice.chems.data <- county_cyp1a1_up %>% 
  group_by(casrn) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(casrn %in% in_chems)

### Load in the in-vitro concentration response data from ICE
ice_data <- get(load("/Volumes/messierkp/Projects/AEP-AOP/LTEA_HepaRG_CYP1A1_up 41 chems for Kyle 220131.RData"))

# split the data by chemical 
ice.data.by.chem <- split(ice_data,as.factor(ice_data$casn))


### Test our own version of fitting the 3-parameter hill model 
# Our goal is to check and compare the values - particularly the parameter 
# standard deviations
source("/Volumes/messierkp/Projects/AEP-AOP/tcpl_my_3hill_fit.R")
source("/Volumes/messierkp/Projects/AEP-AOP/tcpl_3hill_obj.R")

my.3hill.mdl <- my_3hill_fit(ice.data.by.chem[[1]]$logc,ice.data.by.chem[[1]]$resp,log = TRUE)
tcpl.3hill.mdl <- fithill(10^ice.data.by.chem[[1]]$logc,ice.data.by.chem[[1]]$resp)

# This all checks out with tcpl hill function and the curated ICE data 
# And there was much rejoice! 

#### My 2-parameter hill model - testing fit

source("/Volumes/messierkp/Projects/AEP-AOP/tcpl_my_2hill_fit.R")
source("/Volumes/messierkp/Projects/AEP-AOP/tcpl_2hill_obj.R")

my.2hill.mdl <- my_2hill_fit(ice.data.by.chem[[1]]$logc,ice.data.by.chem[[1]]$resp,log = TRUE)

# Make some comparison plots

X <- 10^seq(log10(10^-3),log10(10^3),length.out = 10^3)
logX <- log10(X)

val.3hill <- my.3hill.mdl$par[1] /(1 + 10^( my.3hill.mdl$par[3] * (my.3hill.mdl$par[2] - logX) ))
val.2hill <- my.2hill.mdl$par[1] /(1 + 10^( (my.2hill.mdl$par[2] - logX) ))

df.plot <- data.frame("X" = X, "threeHill" = val.3hill,"twoHill" = val.2hill)

ggplot(df.plot) + geom_line(aes(X,threeHill),color = "red") + 
  geom_line(aes(X,twoHill),color = "blue") + 
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))
