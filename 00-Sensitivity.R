######################################################
# By: Kyle Messier
# Updated, Run, 02/04/2022
# Written in R Version 4.0.2
######################################################
# First script in the sensitivity analysis. This calculates 
# average Css statistics for age, obesity, and httk params to use in their respective
# sensitivity analyses


library(httk)
library(truncnorm)
library(tidyverse)

load_sipes2017()
# load in the saved data
age.by.county <- get(load("/Volumes/SHAG/GeoTox/data/age_by_county_20220201.RData"))
obesity.by.county <- get(load("/Volumes/SHAG/GeoTox/data/obesity_by_county_20220201.RData"))

# css.list <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/httk_css_pre_simulate_20211209.RData"))

# This is the same order of chemicals in the by.county variables such as external and invitro conentration.by,county
in.chems <-  c("98-86-2","92-87-5","92-52-4","117-81-7","133-06-2","532-27-4","133-90-4","57-74-9","510-15-6","94-75-7" ,
               "64-67-5","132-64-9","106-46-7","111-44-4","79-44-7","131-11-3","77-78-1","119-90-4","121-14-2","534-52-1", 
               "51-28-5","121-69-7","107-21-1","51-79-6","76-44-8","822-06-0","77-47-4","123-31-9","72-43-5" , 
               "101-77-9","56-38-2","82-68-8","87-86-5","1120-71-4", "114-26-1","91-22-5","96-09-3","95-80-7","584-84-9" ,
               "95-95-4","1582-09-8")


MC.iter <- 10^3

# Groups 

ages <- c("a","b","c","d","e","f","g","h","i","j","k")
weight <- c("Normal","Obese")

age.cutoffs <- rbind(c(0,2),
                     c(3,5),
                     c(6,10),
                     c(11,15),
                     c(16,20),
                     c(21,30),
                     c(31,40),
                     c(41,50),
                     c(51,60),
                     c(61,70),
                     c(71,100))

group = matrix(NA,nrow = 11*2,ncol = 2)
l=1
for (i in ages){
  for (j in weight){
    
    group[l,] <- c(i,j)
    l=l+1
  }
}  


if (file.exists(file = "/Volumes/SHAG/GeoTox/data/css_list_sensitivity_simulate_20220201.RData")){
  css.list <- get(load("/Volumes/SHAG/GeoTox/data/css_list_sensitivity_simulate_20220201.RData"))
}else{
  css.list <- list()
  for (i in 1:length(in.chems)){
    
    val <- matrix(NA,MC.iter,nrow(group))
    for (j in 1:nrow(group)){
      
      print(c(i,j))
      if (group[j,1]=="a"){
        agelim <- age.cutoffs[1,]
      }else if(group[j,1]=="b"){
        agelim <- age.cutoffs[2,]
      }else if(group[j,1]=="c"){
        agelim <- age.cutoffs[3,]
      }else if(group[j,1]=="d"){
        agelim <- age.cutoffs[4,]
      }else if(group[j,1]=="e"){
        agelim <- age.cutoffs[5,]
      }else if(group[j,1]=="f"){
        agelim <- age.cutoffs[6,]
      }else if(group[j,1]=="g"){
        agelim <- age.cutoffs[7,]
      }else if(group[j,1]=="h"){
        agelim <- age.cutoffs[8,]
      }else if(group[j,1]=="i"){
        agelim <- age.cutoffs[9,]
      }else if(group[j,1]=="j"){
        agelim <- age.cutoffs[10,]
      }else if(group[j,1]=="k"){
        agelim <- age.cutoffs[11,]
      }
      
      httkpoplist = list( method = "vi", 
                          gendernum = NULL,
                          agelim_years = agelim,
                          agelim_months = NULL, 
                          weight_category = group[j,2],
                          reths =c(
                            "Mexican American",
                            "Other Hispanic",
                            "Non-Hispanic White",
                            "Non-Hispanic Black",
                            "Other"
                          ))
      
      mcs <- create_mc_samples(chem.cas = in.chems[i],samples = MC.iter,httkpop.generate.arg.list = httkpoplist)
      css <- calc_analytic_css(chem.cas = in.chems[i],parameters = mcs, model = "3compartmentss")
      
      val[,j] <- css
      
      
      
    }
    
    css.list[[i]] <- val
  }
  
  save(css.list,file = "/Volumes/SHAG/GeoTox/data/css_list_sensitivity_simulate_20220201.RData")
}




css.sensitivity.age <- NULL
css.sensitivity.obesity <- NULL
css.sensitivity.httk <- NULL

for (i in 1:length(age.by.county)){

  css.sensitivity.age[[i]]<- matrix(NA,nrow = MC.iter,ncol = 41)
  css.sensitivity.obesity[[i]]<- matrix(NA,nrow = MC.iter,ncol = 41)
  css.sensitivity.httk[[i]]<- matrix(NA,nrow = MC.iter,ncol = 41)

  idx.age.group1 <- round(age.by.county[[i]])<=2
  idx.age.group2 <- round(age.by.county[[i]])>=3 & round(age.by.county[[i]])<=5
  idx.age.group3 <- round(age.by.county[[i]])>=6 & round(age.by.county[[i]])<=10
  idx.age.group4 <- round(age.by.county[[i]])>=11 & round(age.by.county[[i]])<=15
  idx.age.group5 <- round(age.by.county[[i]])>=16 & round(age.by.county[[i]])<=20
  idx.age.group6 <- round(age.by.county[[i]])>=21 & round(age.by.county[[i]])<=30
  idx.age.group7 <- round(age.by.county[[i]])>=31 & round(age.by.county[[i]])<=40
  idx.age.group8 <- round(age.by.county[[i]])>=41 & round(age.by.county[[i]])<=50
  idx.age.group9 <- round(age.by.county[[i]])>=51 & round(age.by.county[[i]])<=60
  idx.age.group10 <- round(age.by.county[[i]])>=61 & round(age.by.county[[i]])<=70
  idx.age.group11 <- round(age.by.county[[i]])>=71
  
  
  idx.obesity1 <- obesity.by.county[[i]] == "Normal"
  idx.obesity2 <- obesity.by.county[[i]] == "Obese"  
  
  

    seq1 <- seq(1,22,by=2)
    seq2 <- seq(2,22,by=2)
    tot_seq <- rbind(seq1,seq2)
    for (j in 1:41){ # For each chemical
      
      print(c(i,j))
      
      css.medians.by.age <- matrix(NA,nrow = 11)
      for (k in 1:11){
        css.medians.by.age[k] <- median(css.list[[j]][,seq1[k]:seq2[k]])
      }
      
      css.medians.by.obesity <- matrix(NA,nrow = 2)
      for (k in 1:2){
        css.medians.by.obesity[k] <- median(css.list[[j]][,tot_seq[k,]])
      }
      
      
      # Age
      css.sensitivity.age[[i]][idx.age.group1,j] <- css.medians.by.age[1]
      css.sensitivity.age[[i]][idx.age.group2,j] <- css.medians.by.age[2]
      css.sensitivity.age[[i]][idx.age.group3,j] <- css.medians.by.age[3]
      css.sensitivity.age[[i]][idx.age.group4,j] <- css.medians.by.age[4]
      css.sensitivity.age[[i]][idx.age.group5,j] <- css.medians.by.age[5]
      css.sensitivity.age[[i]][idx.age.group6,j] <- css.medians.by.age[6]
      css.sensitivity.age[[i]][idx.age.group7,j] <- css.medians.by.age[7]
      css.sensitivity.age[[i]][idx.age.group8,j] <- css.medians.by.age[8]
      css.sensitivity.age[[i]][idx.age.group9,j] <- css.medians.by.age[9]
      css.sensitivity.age[[i]][idx.age.group10,j] <- css.medians.by.age[10]
      css.sensitivity.age[[i]][idx.age.group11,j] <- css.medians.by.age[11]
      
      # Obesity
      css.sensitivity.obesity[[i]][idx.obesity1,j] <- css.medians.by.obesity[1]
      css.sensitivity.obesity[[i]][idx.obesity2,j] <- css.medians.by.obesity[2]
      
      # HTTK parameters
      median.age <- round(median(age.by.county[[i]]))
      age.category <- sapply(1:nrow(age.cutoffs),function(x) median.age >= age.cutoffs[x,1] & median.age<=age.cutoffs[x,2])
      css.sensitivity.httk[[i]][,j] <- sample(css.list[[j]][,seq1[age.category]],MC.iter,replace = TRUE)

    }


  
}

save(css.sensitivity.obesity,
     file = "/Volumes/SHAG/GeoTox/data/css_by_county_sensitivity_obesity_20220201.RData")

save(css.sensitivity.age,
     file = "/Volumes/SHAG/GeoTox/data/css_by_county_sensitivity_age_20220201.RData")

save(css.sensitivity.httk,
     file = "/Volumes/SHAG/GeoTox/data/css_by_county_sensitivity_httk_20220201.RData")
