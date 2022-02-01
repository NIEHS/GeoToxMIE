census.age.sim <- function(MC.iter, census.age){
#  Create a simulated age distribution by county
#  A utility function for the GeoToxMIE Monte Carlo Analysis
#
  
  
# The census.age data is by 5 year age groups and by county from the census
  age.seq <- 0:89

  uFIPS <- unique(census.age$FIPS)
  uAGE <- unique(census.age$AGEGRP)
  
  age.mat <- matrix(NA,nrow = length(uFIPS),ncol = length(age.seq))

  grp <- seq(0,89,by = 5)  
  for (i in 1:length(uFIPS)){
    age.i <- subset(census.age,FIPS == uFIPS[i])
    total.pop <- subset(age.i,AGEGRP == 0)$TOT_POP
    age.i <- age.i[-1,]
    age.i$prop <- age.i$TOT_POP / total.pop
    
    for (j in age.seq){
      
      # Within each age group, we assume each year is equally probable
      # print(c(i,j))
      idx <- max(age.i$AGEGRP[grp <= j])
      prop.j <- age.i$prop[idx] / 5
      # j+1 because we iterate starting at 0
      age.mat[i,j+1] <- prop.j
    }
  }
  

  sample.fun <- function(x){
    sample(age.seq,size = MC.iter,prob = age.mat[x,],replace = TRUE)
  }
  y <- lapply(1:length(uFIPS), sample.fun)
  
}
