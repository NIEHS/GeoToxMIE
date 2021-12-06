census.age.sim <- function(MC.iter, census.age){
# 
#
#
  

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
      
      # print(c(i,j))
      idx <- max(age.i$AGEGRP[grp <= j])
      prop.j <- age.i$prop[idx] / 5
      age.mat[i,j+1] <- prop.j
    }
  }
  

  sample.fun <- function(x){
    sample(age.seq,size = MC.iter,prob = age.mat[x,],replace = TRUE)
  }
  y <- lapply(1:length(uFIPS), sample.fun)
  
}
