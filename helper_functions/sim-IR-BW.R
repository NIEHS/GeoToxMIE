sim.IR.BW <- function(MC.iter,age.data){
  # Simulate the inhalation rate and body weight based on the age
  # Written: KPM
  # QC, 12/8/21, KPM
  
  
  # Data comes from https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter06.pdf
  # Table 6.7 Distribution percentiles of physiological daily inhalation rates per unit
  # body weight (m3/kg-day) for free living normal weight males and females aged 2 months to 96 years
  age.lower.bound <- c(0,1,2,5,7,11,23,30,40,65)
  age.upper.bound <- c(1,2,5,7,11,23,30,40,65,100)
  
  mean.male <- c(0.495,0.48,0.44,0.42,0.37,0.3,0.25,0.24,0.23,0.19)
  sd.male <- c(0.08,0.06,0.04,0.05,0.06,0.05,0.04,0.03,0.04,0.03)
  
  mean.female <- c(0.48,0.45,0.44,0.40,0.35,0.27,0.23,0.24,0.21,0.17)
  sd.female <- c(0.075,0.08,0.07,0.05,0.06,0.05,0.04,0.04,0.04,0.04)
  
  m1 <- apply(cbind(mean.male,mean.female),1,mean)
  sd1 <- apply(cbind(sd.male,sd.female),1,mean)

  # sim.age.mean.male <- sim.age.mean.female <- sim.age.sd.male <- sim.age.sd.female <- age.data * NA
 
  idx.fun <- function(x){
    IR.mean <- IR.sd <- rep(NA,MC.iter)
    idx1 <- age.data[[x]] >= age.lower.bound[1] & age.data[[x]] < age.upper.bound[1]
    idx2 <- age.data[[x]] >= age.lower.bound[2] & age.data[[x]] < age.upper.bound[2]
    idx3 <- age.data[[x]] >= age.lower.bound[3] & age.data[[x]] < age.upper.bound[3]
    idx4 <- age.data[[x]] >= age.lower.bound[4] & age.data[[x]] < age.upper.bound[4]
    idx5 <- age.data[[x]] >= age.lower.bound[5] & age.data[[x]] < age.upper.bound[5]
    idx6 <- age.data[[x]] >= age.lower.bound[6] & age.data[[x]] < age.upper.bound[6]
    idx7 <- age.data[[x]] >= age.lower.bound[7] & age.data[[x]] < age.upper.bound[7]
    idx8 <- age.data[[x]] >= age.lower.bound[8] & age.data[[x]] < age.upper.bound[8]
    idx9 <- age.data[[x]] >= age.lower.bound[9] & age.data[[x]] < age.upper.bound[9]
    idx10 <- age.data[[x]] >= age.lower.bound[10] & age.data[[x]] < age.upper.bound[10]
    
    IR.mean[idx1] <- m1[1]
    IR.mean[idx2] <- m1[2]
    IR.mean[idx3] <- m1[3]
    IR.mean[idx4] <- m1[4]
    IR.mean[idx5] <- m1[5]
    IR.mean[idx6] <- m1[6]
    IR.mean[idx7] <- m1[7]
    IR.mean[idx8] <- m1[8]
    IR.mean[idx9] <- m1[9]
    IR.mean[idx10] <- m1[10]
    
    IR.sd[idx1] <- sd1[1]
    IR.sd[idx2] <- sd1[2]
    IR.sd[idx3] <- sd1[3]
    IR.sd[idx4] <- sd1[4]
    IR.sd[idx5] <- sd1[5]
    IR.sd[idx6] <- sd1[6]
    IR.sd[idx7] <- sd1[7]
    IR.sd[idx8] <- sd1[8]
    IR.sd[idx9] <- sd1[9]
    IR.sd[idx10] <- sd1[10]
    
    df <- data.frame("IR.mean" = IR.mean,"IR.sd" = IR.sd)
    return(df)
  }
  

  IR <- lapply(1:length(age.data),idx.fun)
  

  IR.fun <- function(x){
    rtruncnorm(length(IR[[x]]$IR.mean),0,Inf,mean = IR[[x]]$IR.mean,sd = IR[[x]]$IR.sd)
  }
  
  IR.by.county <- lapply(1:length(IR),IR.fun)
  
}
