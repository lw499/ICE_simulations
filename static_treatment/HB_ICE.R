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
#ptm <- proc.time()
#subset data
yfitog = glm(Y~ (A+L1), family = binomial(), data = dffull) #This is from the data generation mechanism
paramtmp = (yfitog)$coef
treats = c(0,1)

##################
######time 5######
##################
meany5 = NULL
for(i in 1:2)
{ 
  abar = treats[i]
  y5dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0,]; y5dat$A_4 = abar;
  y5dat$A = y5dat$A_4; y5dat$L1 = y5dat$L1_4;
  y5dat$y5pred = predict(yfitog, newdata = y5dat, type="response"); 
  
  tmp = y5dat[y5dat$Y_4==0 & y5dat$C_4==0,] #this covers those who are alive and who have not been censored
  y5dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0  & tmpdata$Y_3==0,];
  y5fitextra = glm(y5pred ~ L1_3 + A_3 + L1_3*A_3, family = binomial(), data = tmp) ; y5dat$A_3 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y5dat)), y5dat$A_3, y5dat$L1_3)  %*% matrix(paramtmp, nrow=3))
  y5dat$y5pred = predict(y5fitextra, newdata = y5dat, type="response")*(1-predicttmp) + predicttmp
  
  tmp = y5dat[y5dat$Y_3==0 & y5dat$C_3==0,] #this covers those who are alive and who have not been censored
  y5dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0  & tmpdata$Y_2==0,];
  y5fitextra = glm(y5pred ~ L1_2 + A_2 + L1_2*A_2, family = binomial(), data = tmp) ; y5dat$A_2 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y5dat)), y5dat$A_2, y5dat$L1_2)  %*% matrix(paramtmp, nrow=3))
  y5dat$y5pred = predict(y5fitextra, newdata = y5dat, type="response")*(1-predicttmp) + predicttmp
  
  tmp = y5dat[y5dat$Y_2==0 & y5dat$C_2==0,] #this covers those who are alive and who have not been censored
  y5dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0  & tmpdata$Y_1==0,];
  y5fitextra = glm(y5pred ~ L1_1 + A_1 + L1_1*A_1, family = binomial(), data = tmp) ; y5dat$A_1 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y5dat)), y5dat$A_1, y5dat$L1_1)  %*% matrix(paramtmp, nrow=3))
  y5dat$y5pred = predict(y5fitextra, newdata = y5dat, type="response")*(1-predicttmp) + predicttmp
  
  tmp = y5dat[y5dat$Y_1==0 & y5dat$C_1==0,] #this covers those who are alive and who have not been censored
  y5dat = tmpdata
  y5fitextra = glm(y5pred ~ L1_0 + A_0 + L1_0*A_0, family = binomial(), data = tmp) ; y5dat$A_0 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y5dat)), y5dat$A_0, y5dat$L1_0)  %*% matrix(paramtmp, nrow=3))
  y5dat$y5pred = predict(y5fitextra, newdata = y5dat, type="response")*(1-predicttmp) + predicttmp
  
  meany5tmp = c(i,mean(y5dat$y5pred))
  
  meany5 = rbind(meany5, meany5tmp)
}
#proc.time() - ptm

##################
######time 4######
##################
meany4 = NULL
for(i in 1:2)
{ 
  abar = treats[i]
  y4dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0,]; y4dat$A_3 = abar;
  y4dat$A = y4dat$A_3; y4dat$L1 = y4dat$L1_3;
  y4dat$y4pred = predict(yfitog, newdata = y4dat, type="response"); 
  
  tmp = y4dat[y4dat$Y_3==0 & y4dat$C_3==0,] #this covers those who are alive and who have not been censored
  y4dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0  & tmpdata$Y_2==0,];
  y4fitextra = glm(y4pred ~ L1_2 + A_2 + L1_2*A_2, family = binomial(), data = tmp) ; y4dat$A_2 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y4dat)), y4dat$A_2, y4dat$L1_2)  %*% matrix(paramtmp, nrow=3))
  y4dat$y4pred = predict(y4fitextra, newdata = y4dat, type="response")*(1-predicttmp) + predicttmp
  
  tmp = y4dat[y4dat$Y_2==0 & y4dat$C_2==0,] #this covers those who are alive and who have not been censored
  y4dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0  & tmpdata$Y_1==0,];
  y4fitextra = glm(y4pred ~ L1_1 + A_1 + L1_1*A_1, family = binomial(), data = tmp) ; y4dat$A_1 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y4dat)), y4dat$A_1, y4dat$L1_1)  %*% matrix(paramtmp, nrow=3))
  y4dat$y4pred = predict(y4fitextra, newdata = y4dat, type="response")*(1-predicttmp) + predicttmp
  
  tmp = y4dat[y4dat$Y_1==0 & y4dat$C_1==0,] #this covers those who are alive and who have not been censored
  y4dat = tmpdata
  y4fitextra = glm(y4pred ~ L1_0 + A_0 + L1_0*A_0, family = binomial(), data = tmp) ; y4dat$A_0 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y4dat)), y4dat$A_0, y4dat$L1_0)  %*% matrix(paramtmp, nrow=3))
  y4dat$y4pred = predict(y4fitextra, newdata = y4dat, type="response")*(1-predicttmp) + predicttmp
  
  meany4tmp = c(i,mean(y4dat$y4pred))
  
  meany4 = rbind(meany4, meany4tmp)
}

##################
######time 3######
##################
meany3 = NULL
for(i in 1:2)
{ 
  abar = treats[i]
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0,]; y3dat$A_2 = abar;
  y3dat$A = y3dat$A_2; y3dat$L1 = y3dat$L1_2;
  y3dat$y3pred = predict(yfitog, newdata = y3dat, type="response"); 
  
  tmp = y3dat[y3dat$Y_2==0 & y3dat$C_2==0,] #this covers those who are alive and who have not been censored
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0  & tmpdata$Y_1==0,];
  y3fitextra = glm(y3pred ~ L1_1 + A_1 + L1_1*A_1, family = binomial(), data = tmp) ; y3dat$A_1 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y3dat)), y3dat$A_1, y3dat$L1_1)  %*% matrix(paramtmp, nrow=3))
  y3dat$y3pred = predict(y3fitextra, newdata = y3dat, type="response")*(1-predicttmp) + predicttmp
  
  tmp = y3dat[y3dat$Y_1==0 & y3dat$C_1==0,] #this covers those who are alive and who have not been censored
  y3dat = tmpdata
  y3fitextra = glm(y3pred ~ L1_0 + A_0 + L1_0*A_0, family = binomial(), data = tmp) ; y3dat$A_0 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y3dat)), y3dat$A_0, y3dat$L1_0)  %*% matrix(paramtmp, nrow=3))
  y3dat$y3pred = predict(y3fitextra, newdata = y3dat, type="response")*(1-predicttmp) + predicttmp
  
  meany3tmp = c(i,mean(y3dat$y3pred))
  
  meany3 = rbind(meany3, meany3tmp)
}

##################
######time 2######
##################
meany2 = NULL
for(i in 1:2)
{ 
  abar = treats[i]
  y2dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0,]; y2dat$A_1 = abar;
  y2dat$A = y2dat$A_1; y2dat$L1 = y2dat$L1_1;
  y2dat$y2pred = predict(yfitog, newdata = y2dat, type="response"); 
 
  tmp = y2dat[y2dat$Y_1==0 & y2dat$C_1==0,] #this covers those who are alive and who have not been censored
  y2dat = tmpdata
  y2fitextra = glm(y2pred ~ L1_0 + A_0 + L1_0*A_0, family = binomial(), data = tmp) ; y2dat$A_0 = abar;
      predicttmp = plogis(cbind(rep(1, nrow(y2dat)), y2dat$A_0, y2dat$L1_0)  %*% matrix(paramtmp, nrow=3))
  y2dat$y2pred = predict(y2fitextra, newdata = y2dat, type="response")*(1-predicttmp) + predicttmp
  
  meany2tmp = c(i,mean(y2dat$y2pred))
  
  meany2 = rbind(meany2, meany2tmp)
}
  
##################
######time 1######
##################
meany1 = NULL
for(i in 1:2)
{ 
  abar = treats[i]
  y1dat = tmpdata; y1dat$A_0 = abar;
  y1dat$A = y1dat$A_0; y1dat$L1 = y1dat$L1_0;
  y1dat$y1pred = predict(yfitog, newdata = y1dat, type="response"); 
  
  meany1tmp = c(i,mean(y1dat$y1pred))
  
  meany1 = rbind(meany1, meany1tmp)
}

myparam = cbind(meany1, meany2, meany3, meany4, meany5)
#proc.time() - ptm

return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"HB_ICE.csv")

stopCluster(cl)