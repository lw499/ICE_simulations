library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(detectCores()-2)                                              
registerDoParallel(cl)                                                                          
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  options(warn=-1)
  library(geepack);library(MASS);library(ResourceSelection);library(ltmle); library(SuperLearner)
  library(dplyr)
  library(data.table)
  #library(reshape2)  #do not use for data frame only
  
  logit <- function(term) {
    return( ifelse(!is.na(term),log(term/(1-term)),NA) )
  }
  
  EXPIT <- function(term) {
    return( ifelse(!is.na(term),exp(term)/(1+exp(term)),NA) )
  }
  
  source("datagen.R")
  set.seed(1129)
  seeds = floor(runif(1000)*10^8);
  set.seed(seeds[m])
  
  n <- 1000
  K <- 5
  alpha0=-1.5; alpha1=1; alpha2=0; alpha3=0;
  beta0=-2; beta1=-2; beta2=0; beta3=0;
  theta0=-2; theta1=-2; theta2=1; theta3=0; theta4=0;
  eta0=-3; eta1=-1; eta2=0.75; eta3=0; eta4=0;
  gamma0=-2; gamma1=-1; gamma2=-1; gamma3=0; gamma4=0
  sigma=0.1
  df <- lapply(as.list(1:n), FUN=function(ind){
    datagen(ind, K=K,
            alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3,
            beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3,
            theta0=theta0, theta1=theta1, theta2=theta2, theta3=theta3, theta4=theta4,
            eta0=eta0, eta1=eta1, eta2=eta2, eta3=eta3, eta4=eta4, sigma=sigma)
  })
  
  dffull <- rbindlist(df)
  
  
  dffull[, paste("lag_A") := shift(A, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_C") := shift(C, 1, NA, type='lag'), by=id]
  
  
  afit = glm(A ~ L1, family = binomial(), data = dffull[!is.na(dffull$lag_A) & dffull$lag_A==0,]) #This is from the data generation mechanism
  cfit = glm(C ~ A + L1, family = binomial(), data = dffull) #This is from the data generation mechanism
  
  dffull$pred_obs = predict(afit, newdata = dffull, type="response")
  dffull$pred_obs = ifelse(dffull$A==1, dffull$pred_obs, 1-dffull$pred_obs)
    dffull$pred_obs = ifelse(!is.na(dffull$lag_A) & dffull$lag_A==1, 1, dffull$pred_obs)
  dffull$pred_obsc = predict(cfit, newdata = dffull, type="response")
  dffull$pred_obsc = ifelse(dffull$C==1, NA, 1-dffull$pred_obsc)
  
  #system.time( dcast(dffull, id ~ t0, value.var = c("L1","L2","A","C","Y")) )
  dffullwide = dcast(dffull, id ~ t0, value.var = c("L1","L2","A","C","Y","pred_obs","pred_obsc"))
  
  tmpdata = dffullwide
  tmpdata$C_1 = ifelse(is.na(tmpdata$Y_0) & tmpdata$C_0==1,1,tmpdata$C_1)
  tmpdata$C_2 = ifelse(is.na(tmpdata$Y_1) & tmpdata$C_1==1,1,tmpdata$C_2)
  tmpdata$C_3 = ifelse(is.na(tmpdata$Y_2) & tmpdata$C_2==1,1,tmpdata$C_3)
  tmpdata$C_4 = ifelse(is.na(tmpdata$Y_3) & tmpdata$C_3==1,1,tmpdata$C_4)
  
  tmpdata$Y_1 = ifelse(tmpdata$Y_0==1,1,tmpdata$Y_1)
  tmpdata$Y_2 = ifelse(!is.na(tmpdata$Y_1) & tmpdata$Y_1==1,1,tmpdata$Y_2)
  tmpdata$Y_3 = ifelse(!is.na(tmpdata$Y_2) & tmpdata$Y_2==1,1,tmpdata$Y_3)
  tmpdata$Y_4 = ifelse(!is.na(tmpdata$Y_3) & tmpdata$Y_3==1,1,tmpdata$Y_4)
  
  tmpdata$Y_5 = tmpdata$Y_4
  tmpdata$Y_4 = tmpdata$Y_3
  tmpdata$Y_3 = tmpdata$Y_2
  tmpdata$Y_2 = tmpdata$Y_1
  tmpdata$Y_1 = tmpdata$Y_0
  
  tmpdata$C_5 = tmpdata$C_4
  tmpdata$C_4 = tmpdata$C_3
  tmpdata$C_3 = tmpdata$C_2
  tmpdata$C_2 = tmpdata$C_1
  tmpdata$C_1 = tmpdata$C_0
  
  tmpdata$Y_0 = tmpdata$C_0 = NULL;
  
  tmpdata$id = seq(1,n,by=1)
  tmpdata$pi6 <- tmpdata$pi5 <- tmpdata$pi4 <- tmpdata$pi3 <- tmpdata$pi2 <- tmpdata$pi1 <- tmpdata$pi0 <- NA
  tmpdata$pi6c <- tmpdata$pi5c <- tmpdata$pi4c <- tmpdata$pi3c <- tmpdata$pi2c <- tmpdata$pi1c <- tmpdata$pi0c <- NA
  
  tmpdata$pi0 = tmpdata$pred_obs_0
  tmpdata$pi1 = tmpdata$pi0*tmpdata$pred_obs_1
  tmpdata$pi2 = tmpdata$pi1*tmpdata$pred_obs_2
  tmpdata$pi3 = tmpdata$pi2*tmpdata$pred_obs_3
  tmpdata$pi4 = tmpdata$pi3*tmpdata$pred_obs_4
#  tmpdata$pi5 = tmpdata$pi4*tmpdata$pred_obs_5
#  tmpdata$pi6 = tmpdata$pi5*tmpdata$pred_obs_6
  
  tmpdata$pi0c = tmpdata$pred_obsc_0
  tmpdata$pi1c = tmpdata$pi0c*tmpdata$pred_obsc_1
  tmpdata$pi2c = tmpdata$pi1c*tmpdata$pred_obsc_2
  tmpdata$pi3c = tmpdata$pi2c*tmpdata$pred_obsc_3
  tmpdata$pi4c = tmpdata$pi3c*tmpdata$pred_obsc_4
#  tmpdata$pi5c = tmpdata$pi4c*tmpdata$pred_obsc_5
#  tmpdata$pi6c = tmpdata$pi5c*tmpdata$pred_obsc_6
  
  ##calculate risk
  mean = NULL
    ind = NA;
  a=1
  #time 1
  tmp = tmpdata[tmpdata$A_0==a & tmpdata$C_1==0,]
  component1 = tmp$Y_1/tmp$pi0; component1c = tmp$Y_1/tmp$pi0c;
  comp1 = sum(component1*component1c); 
  if(nrow(tmp)>0){param1 = comp1/n} else{param1=NA}
  #time 2
  tmp  = tmpdata[tmpdata$A_0==a &tmpdata$A_1==a & tmpdata$Y_1==0 & tmpdata$C_2==0,]
  component2 = tmp$Y_2/tmp$pi1; component2c = tmp$Y_2/tmp$pi1c;
  comp2 = sum(component2*component2c)
  if(nrow(tmp)>0){param2 = (comp1 + comp2)/n} else{param2=NA}
  #time 3
  tmp  = tmpdata[tmpdata$A_0==a & tmpdata$A_1==a & tmpdata$A_2==a & tmpdata$Y_2==0 & tmpdata$C_3==0,]
  component3 = tmp$Y_3/tmp$pi2; component3c = tmp$Y_3/tmp$pi2c;
  comp3 = sum(component3*component3c)
  if(nrow(tmp)>0){param3 = (comp1 + comp2 + comp3)/n} else{param3=NA}
  #time 4
  tmp  = tmpdata[tmpdata$A_0==a & tmpdata$A_1==a & tmpdata$A_2==a & tmpdata$A_3==a & tmpdata$Y_3==0 & tmpdata$C_4==0,]
  component4 = tmp$Y_4/tmp$pi3; component4c = tmp$Y_4/tmp$pi3c;
  comp4 = sum(component4*component4c)
  if(nrow(tmp)>0){param4 = (comp1 + comp2 + comp3 + comp4)/n} else{param4=NA}
  #time 5
  tmp  = tmpdata[tmpdata$A_0==a & tmpdata$A_1==a & tmpdata$A_2==a & tmpdata$A_3==a & tmpdata$A_4==a & tmpdata$Y_4==0 & tmpdata$C_5==0,]
  component5 = tmp$Y_5/tmp$pi4; component5c = tmp$Y_5/tmp$pi4c;
  comp5 = sum(component5*component5c)
  if(nrow(tmp)>0){param5 = (comp1 + comp2 + comp3 + comp4 + comp5)/n} else{param5=NA}

  ind = ifelse(comp1==0 | comp2==0 | comp3==0 | comp4==0 | comp5==0, 1, ind)
  
  mean = c(param1, param2, param3, param4, param5, ind)
  myparam = mean
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"ipw.csv")

stopCluster(cl)
