######################################################
# By: Kyle Messier, Kristin Eccles
# Date: Feb 1, 2022
# Written in R Version 4.0.2
######################################################
# Fourth script in the main CYCP1A1 analysis pipeline
# This script takes the precalculated css values and does
# random sampling by county to produce a valid monte carlo iteration
# of Css by county



library(httk)
library(truncnorm)
library(tidyverse)

# load in the saved data
age.by.county <- get(load("/Volumes/SHAG/GeoTox/data/age_by_county_20220201.RData"))
obesity.by.county <- get(load("/Volumes/SHAG/GeoTox/data/obesity_by_county_20220201.RData"))

css.list <- get(load("/Volumes/SHAG/GeoTox/data/httk_css_pre_simulate_20220201.RData"))


MC.iter <- 10^3

# Groups 

ages <- c("a","b","c","d","e","f","g","h","i","j","k")
weight <- c("Normal","Obese")


group = matrix(NA,nrow = 11*2,ncol = 2)
l=1
for (i in ages){
  for (j in weight){
    
    group[l,] <- c(i,j)
    l=l+1
  }
}  

# print(c(i,j))

css.by.county <- NULL
for (i in 1:length(age.by.county)){

  css.by.county[[i]]<- matrix(NA,nrow = MC.iter,ncol = 41)
  
  
  
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
    
    for (j in 1:41){ # For each chemical
      print(c(i,j))
      # Group 1
      idx.1 <-idx.age.group1 & idx.obesity1
      css.by.county[[i]][idx.1,j] <- sample(css.list[[j]][,1],sum(idx.1),replace = TRUE)
      # Group 2
      idx.2 <-idx.age.group1 & idx.obesity2
      css.by.county[[i]][idx.2,j] <- sample(css.list[[j]][,2],sum(idx.2),replace = TRUE)
      # Group 3
      idx.3 <-idx.age.group2 & idx.obesity1
      css.by.county[[i]][idx.3,j] <- sample(css.list[[j]][,3],sum(idx.3),replace = TRUE)
      # Group 4
      idx.4 <-idx.age.group2 & idx.obesity2
      css.by.county[[i]][idx.4,j] <- sample(css.list[[j]][,4],sum(idx.4),replace = TRUE)
      # Group 5
      idx.5 <-idx.age.group3 & idx.obesity1
      css.by.county[[i]][idx.5,j] <- sample(css.list[[j]][,5],sum(idx.5),replace = TRUE)
      # Group 6
      idx.6 <-idx.age.group3 & idx.obesity2
      css.by.county[[i]][idx.6,j] <- sample(css.list[[j]][,6],sum(idx.6),replace = TRUE)
      # Group 7
      idx.7 <-idx.age.group4 & idx.obesity1
      css.by.county[[i]][idx.7,j] <- sample(css.list[[j]][,7],sum(idx.7),replace = TRUE)
      # Group 8
      idx.8 <-idx.age.group4 & idx.obesity2
      css.by.county[[i]][idx.8,j] <- sample(css.list[[j]][,8],sum(idx.8),replace = TRUE)
      # Group 9
      idx.9 <-idx.age.group5 & idx.obesity1
      css.by.county[[i]][idx.9,j] <- sample(css.list[[j]][,9],sum(idx.9),replace = TRUE)
      # Group 10
      idx.10 <-idx.age.group5 & idx.obesity2
      css.by.county[[i]][idx.10,j] <- sample(css.list[[j]][,10],sum(idx.10),replace = TRUE)
      # Group 11
      idx.11 <-idx.age.group6 & idx.obesity1
      css.by.county[[i]][idx.11,j] <- sample(css.list[[j]][,11],sum(idx.11),replace = TRUE)
      # Group 12
      idx.12 <-idx.age.group6 & idx.obesity2
      css.by.county[[i]][idx.12,j] <- sample(css.list[[j]][,12],sum(idx.12),replace = TRUE)
      # Group 13
      idx.13 <-idx.age.group7 & idx.obesity1
      css.by.county[[i]][idx.13,j] <- sample(css.list[[j]][,13],sum(idx.13),replace = TRUE)
      # Group 14
      idx.14 <-idx.age.group7 & idx.obesity2
      css.by.county[[i]][idx.14,j] <- sample(css.list[[j]][,14],sum(idx.14),replace = TRUE)
      # Group 15
      idx.15 <-idx.age.group8 & idx.obesity1
      css.by.county[[i]][idx.15,j] <- sample(css.list[[j]][,15],sum(idx.15),replace = TRUE)
      # Group 16
      idx.16 <-idx.age.group8 & idx.obesity2
      css.by.county[[i]][idx.16,j] <- sample(css.list[[j]][,16],sum(idx.16),replace = TRUE)
      # Group 17
      idx.17 <-idx.age.group9 & idx.obesity1
      css.by.county[[i]][idx.17,j] <- sample(css.list[[j]][,17],sum(idx.17),replace = TRUE)
      # Group 18
      idx.18 <-idx.age.group9 & idx.obesity2
      css.by.county[[i]][idx.18,j] <- sample(css.list[[j]][,18],sum(idx.18),replace = TRUE)
      # Group 19
      idx.19 <-idx.age.group10 & idx.obesity1
      css.by.county[[i]][idx.19,j] <- sample(css.list[[j]][,19],sum(idx.19),replace = TRUE)
      # Group 20
      idx.20 <-idx.age.group10 & idx.obesity2
      css.by.county[[i]][idx.20,j] <- sample(css.list[[j]][,20],sum(idx.20),replace = TRUE)
      # Group 21
      idx.21 <-idx.age.group11 & idx.obesity1
      css.by.county[[i]][idx.21,j] <- sample(css.list[[j]][,21],sum(idx.21),replace = TRUE)
      # Group 22
      idx.22 <-idx.age.group11 & idx.obesity2
      css.by.county[[i]][idx.22,j] <- sample(css.list[[j]][,22],sum(idx.22),replace = TRUE)
      
    }


  
}

save(css.by.county,file = "/Volumes/SHAG/GeoTox/data/css_by_county_20220201.RData")
