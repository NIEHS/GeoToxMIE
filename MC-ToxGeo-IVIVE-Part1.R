# Chemicals in CYP1A1 and NATA
in.chems <-  c("98-86-2","92-87-5","92-52-4","117-81-7","133-06-2","532-27-4","133-90-4","57-74-9","510-15-6","94-75-7" ,
 "64-67-5","132-64-9","106-46-7","111-44-4","79-44-7","131-11-3","77-78-1","119-90-4","121-14-2","534-52-1", 
 "51-28-5","121-69-7","107-21-1","51-79-6","76-44-8","822-06-0","77-47-4","123-31-9","72-43-5" , 
 "101-77-9","56-38-2","82-68-8","87-86-5","1120-71-4", "114-26-1","91-22-5","96-09-3","95-80-7","584-84-9" ,
 "95-95-4","1582-09-8")

# Create the groups where each row represents a unique
# combination of age, weight, and kidney categories
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

# Preallocate a list
css.list <- vector(mode = "list", length = length(in.chems))

# Do the loops - by chemical, by group
MC.iter <- 1000
val <- matrix(NA,nrow = MC.iter,ncol = 10)

for (i in 1:length(in.chems)){
  
  for (j in 1:nrow(group)){
    
    print(c(i,j))
    if (group[j,1]=="a"){
      agelim <- c(0,5)
    }else if(group[j,1]=="b"){
      agelim <- c(6,11)
    }else if(group[j,1]=="c"){
      agelim <- c(12,19)
    }else if(group[j,1]=="d"){
      agelim <- c(20,65)
    }else if(group[j,1]=="e"){
      agelim <- c(66,79)
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


save(css.list,file = "/Volumes/messierkp/Projects/AEP-AOP/httk_css_pre_simulate.RData")
  
save(css.list,file = "/Volumes/SHAG/GeoTox/data/httk_IVIVE/httk_css_pre_simulate.RData")

