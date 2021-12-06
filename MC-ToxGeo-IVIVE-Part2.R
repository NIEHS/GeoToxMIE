


library(httk)
library(truncnorm)
library(tidyverse)

# load in the saved data
age.by.county <- get(load("/Volumes/SHAG/GeoTox/data/age_by_county_20211118.RData"))
obesity.by.county <- get(load("/Volumes/SHAG/GeoTox/data/obesity_by_county_20211118.RData"))

css.list <- get(load("/Volumes/SHAG/GeoTox/data/httk_IVIVE/httk_css_pre_simulate.RData"))


MC.iter <- 10^3

# Groups 

ages <- c("a","b","c","d","e")
weight <- c("Normal","Obese")


group = matrix(NA,nrow = 5*2,ncol = 2)
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
  
    idx.age.group1 <- age.by.county[[i]]<=5
    idx.age.group2 <- age.by.county[[i]]>=6 & age.by.county[[i]]<=11
    idx.age.group3 <- age.by.county[[i]]>=12 & age.by.county[[i]]<=19
    idx.age.group4 <- age.by.county[[i]]>=20 & age.by.county[[i]]<=65
    idx.age.group5 <- age.by.county[[i]]>=66
    
    idx.obesity1 <- obesity.by.county[[i]] == "Normal"
    idx.obesity2 <- obesity.by.county[[i]] == "Obese"
    
    for (j in 1:41){ # For each chemical
      print(c(i,j))
      # Group 1
      idx.1 <-idx.age.group1 & idx.obesity1
      css.by.county[[i]][idx.1,] <- sample(css.list[[j]][,1],sum(idx.1),replace = TRUE)
      # Group 2
      idx.2 <-idx.age.group1 & idx.obesity2
      css.by.county[[i]][idx.2,] <- sample(css.list[[j]][,2],sum(idx.2),replace = TRUE)
      # Group 3
      idx.3 <-idx.age.group2 & idx.obesity1
      css.by.county[[i]][idx.3,] <- sample(css.list[[j]][,3],sum(idx.3),replace = TRUE)
      # Group 4
      idx.4 <-idx.age.group2 & idx.obesity2
      css.by.county[[i]][idx.4,] <- sample(css.list[[j]][,4],sum(idx.4),replace = TRUE)
      # Group 5
      idx.5 <-idx.age.group3 & idx.obesity1
      css.by.county[[i]][idx.5,] <- sample(css.list[[j]][,5],sum(idx.5),replace = TRUE)
      # Group 6
      idx.6 <-idx.age.group3 & idx.obesity2
      css.by.county[[i]][idx.6,] <- sample(css.list[[j]][,6],sum(idx.6),replace = TRUE)
      # Group 7
      idx.7 <-idx.age.group4 & idx.obesity1
      css.by.county[[i]][idx.7,] <- sample(css.list[[j]][,7],sum(idx.7),replace = TRUE)
      # Group 8
      idx.8 <-idx.age.group4 & idx.obesity2
      css.by.county[[i]][idx.8,] <- sample(css.list[[j]][,8],sum(idx.8),replace = TRUE)
      # Group 9
      idx.9 <-idx.age.group5 & idx.obesity1
      css.by.county[[i]][idx.9,] <- sample(css.list[[j]][,9],sum(idx.9),replace = TRUE)
      # Group 10
      idx.10 <-idx.age.group5 & idx.obesity2
      css.by.county[[i]][idx.10,] <- sample(css.list[[j]][,10],sum(idx.10),replace = TRUE)
      
    }


  
}

save(css.by.county,file = "/Volumes/SHAG/GeoTox/data/httk_IVIVE/css_by_county_20211201.RData")
