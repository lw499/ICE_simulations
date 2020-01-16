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
  #system.time( dcast(dffull, id ~ t0, value.var = c("L1","L2","A","C","Y")) )
  dffullwide = dcast(dffull, id ~ t0, value.var = c("L1","L2","A","C","Y"))
  
  # Final simulated dataset
  #system.time(for(i in 1:6)
  #{
  #dffull[, paste("lag.Y",i, sep="") := shift(Y, i, NA, type='lag'), by=id]
  #dffull[, paste("lag.C",i, sep="") := shift(C, i, NA, type='lag'), by=id]
  #dffull[, paste("lag.A",i, sep="") := shift(A, i, NA, type='lag'), by=id]
  #dffull[, paste("lag.L1",i, sep="") :=shift(L1, i, NA,type='lag'),by=id]
  #dffull[, paste("lag.L2",i, sep="") :=shift(L2, i, NA, type='lag'),by=id]
  #})
  
  tmpdata = dffullwide
  #subset data
  yfitog = glm(Y~ (A+L1), family = binomial(), data = dffull) #This is from the data generation mechanism
  paramtmp = (yfitog)$coef
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
  
  y5dat = tmpdata[!is.na(tmpdata$C_5) & tmpdata$C_5==0 & 
                    tmpdata$A_0==1 & tmpdata$A_1==1 &  
                    tmpdata$A_2==1 & tmpdata$A_3==1 &  
                    tmpdata$A_4==1,]; 
  y5fito = glm(Y_5 ~ L1_4, family = binomial(), data = y5dat) ;
  
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==1 & tmpdata$A_1==1 &  
                    tmpdata$A_2==1 &  tmpdata$A_3==1,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==1 & tmpdata$A_1==1 &  
                    tmpdata$A_2==1,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==1 & tmpdata$A_1==1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==1,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  ##################
  ######time 5######
  ##################
  y5dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 & tmpdata$Y_4==0 & 
                     tmpdata$A_0==1 & tmpdata$A_1==1 &  
                     tmpdata$A_2==1 &  tmpdata$A_3==1,];
  y5dat$y5pred = predict(y5fito, newdata = y5dat, type="response")
  
  y5fit = glm(y5pred ~ L1_3, family = binomial(), data = y5dat) ; 
  y5dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & tmpdata$Y_3==0 &
                    tmpdata$A_0==1 & tmpdata$A_1==1 &  
                    tmpdata$A_2==1,];
  y5dat$y5pred = predict(y5fit, newdata = y5dat, type="response")*(1-predict(y4fito, newdata = y5dat, type="response"))+predict(y4fito, newdata = y5dat, type="response"); 

  y5fit = glm(y5pred ~ L1_2, family = binomial(), data = y5dat) ; 
  y5dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==0 &
                    tmpdata$A_0==1 & tmpdata$A_1==1,];
  y5dat$y5pred = predict(y5fit, newdata = y5dat, type="response")*(1-predict(y3fito, newdata = y5dat, type="response"))+predict(y3fito, newdata = y5dat, type="response"); 

  y5fit = glm(y5pred ~ L1_1, family = binomial(), data = y5dat) ; 
  y5dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==0 & tmpdata$A_0==1,]; 
  y5dat$y5pred = predict(y5fit, newdata = y5dat, type="response")*(1-predict(y2fito, newdata = y5dat, type="response"))+predict(y2fito, newdata = y5dat, type="response"); 
 
  y5fit = glm(y5pred ~ L1_0, family = binomial(), data = y5dat) ; 
  y5dat = tmpdata
  y5dat$y5pred = predict(y5fit, newdata = y5dat, type="response")*(1-predict(y1fito, newdata = y5dat, type="response"))+predict(y1fito, newdata = y5dat, type="response"); 
  
  meany5tmp = c(mean(y5dat$y5pred))
  
  meany5 = (meany5tmp)
  
  
  ##################
  ######time 4######
  ##################
  y4dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 &  tmpdata$Y_3==0 & 
                    tmpdata$A_0==1 & tmpdata$A_1==1 &  
                    tmpdata$A_2==1,]; 
  y4dat$y4pred = predict(y4fito, newdata = y4dat, type="response"); 

  y4fit = glm(y4pred ~ L1_2, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==0 &
                    tmpdata$A_0==1 & tmpdata$A_1==1,]; 
  y4dat$y4pred = predict(y4fit, newdata = y4dat, type="response")*(1-predict(y3fito, newdata = y4dat, type="response"))+predict(y3fito, newdata = y4dat, type="response"); 

  y4fit = glm(y4pred ~ L1_1, family = binomial(), data = y4dat) ; 
  y4dat =  tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==0 & tmpdata$A_0==1,]; 
  y4dat$y4pred = predict(y4fit, newdata = y4dat, type="response")*(1-predict(y2fito, newdata = y4dat, type="response"))+predict(y2fito, newdata = y4dat, type="response"); 
  
  y4fit = glm(y4pred ~ L1_0, family = binomial(), data = y4dat) ; 
  y4dat = tmpdata; 
  y4dat$y4pred = predict(y4fit, newdata = y4dat, type="response")*(1-predict(y1fito, newdata = y4dat, type="response"))+predict(y1fito, newdata = y4dat, type="response"); 
  
  meany4tmp = c(mean(y4dat$y4pred))
  
  meany4 = (meany4tmp)
    
  ##################
  ######time 3######
  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==0 & 
                    tmpdata$A_0==1 & tmpdata$A_1==1,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 

  y3fit = glm(y3pred ~ L1_1, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==0 & tmpdata$A_0==1,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*(1-predict(y2fito, newdata = y3dat, type="response"))+predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*(1-predict(y1fito, newdata = y3dat, type="response"))+predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3 = (meany3tmp)
  
  ##################
  ######time 2######
  ##################
  y2dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==0 & tmpdata$A_0==1,]; 
  y2dat$y2pred = predict(y2fito, newdata = y2dat, type="response"); 

  y2fit = glm(y2pred ~ L1_0, family = binomial(), data = y2dat) ; 
  y2dat = tmpdata; 
  y2dat$y2pred = predict(y2fit, newdata = y2dat, type="response")*(1-predict(y1fito, newdata = y2dat, type="response"))+predict(y1fito, newdata = y2dat, type="response"); 
  
  meany2tmp = c(mean(y2dat$y2pred))
  
  meany2 = (meany2tmp)
  
  ##################
  ######time 1######
  ##################
  y1dat = tmpdata; 
  y1dat$y1pred = predict(y1fito, newdata = y1dat, type="response"); 
  
  meany1tmp = c(mean(y1dat$y1pred))
  
  meany1 = (meany1tmp)
  
  
  myparam = cbind(meany1, meany2, meany3, meany4, meany5)
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"HB_ICE_strat.csv")

stopCluster(cl)