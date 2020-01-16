library(doParallel)
library(foreach)

# Calculate the number of cores
getDoParWorkers()                                          
detectCores()                                                      
cl=makeCluster(10)                                              
registerDoParallel(cl)                                                                          
getDoParWorkers()   

#for(m in 1:sim)
myfunc = function(m)
{
  options(warn=-1)
  library(geepack);library(MASS);library(ResourceSelection);library(ltmle); library(SuperLearner)
  library(dplyr)
  library(data.table)
  setDTthreads(1)
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
  dffull[, paste("lag_A") := shift(A, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L") := shift(L1, 1, NA, type='lag'), by=id]
  
  fulldata = dffullwide
  fulldata$C_1 = ifelse(is.na(fulldata$Y_0) & fulldata$C_0==1,1,fulldata$C_1)
  fulldata$C_2 = ifelse(is.na(fulldata$Y_1) & fulldata$C_1==1,1,fulldata$C_2)
  fulldata$C_3 = ifelse(is.na(fulldata$Y_2) & fulldata$C_2==1,1,fulldata$C_3)
  fulldata$C_4 = ifelse(is.na(fulldata$Y_3) & fulldata$C_3==1,1,fulldata$C_4)
  
  fulldata$Y_1 = ifelse(fulldata$Y_0==1,1,fulldata$Y_1)
  fulldata$Y_2 = ifelse(!is.na(fulldata$Y_1) & fulldata$Y_1==1,1,fulldata$Y_2)
  fulldata$Y_3 = ifelse(!is.na(fulldata$Y_2) & fulldata$Y_2==1,1,fulldata$Y_3)
  fulldata$Y_4 = ifelse(!is.na(fulldata$Y_3) & fulldata$Y_3==1,1,fulldata$Y_4)
  
  fulldata$Y_5 = 1-fulldata$Y_4
  fulldata$Y_4 = 1-fulldata$Y_3
  fulldata$Y_3 = 1-fulldata$Y_2
  fulldata$Y_2 = 1-fulldata$Y_1
  fulldata$Y_1 = 1-fulldata$Y_0
  
  fulldata$C_5 = fulldata$C_4
  fulldata$C_4 = fulldata$C_3
  fulldata$C_3 = fulldata$C_2
  fulldata$C_2 = fulldata$C_1
  fulldata$C_1 = fulldata$C_0
  
  fulldata$Y_0 = fulldata$C_0 = NULL;
  
  ## Remove rows that aren't compatible with the treatment regime.
  #fulldata$out = ifelse(fulldata$L1_1==0 & fulldata$A_1==1,1,0)
  #fulldata$out = ifelse(fulldata$L1_0==0 & fulldata$A_0==1,1,fulldata$out)
  #fulldata$out = ifelse(fulldata$L1_0==1 & fulldata$A_1==0,1,fulldata$out)
  
  #tmpdata2 = fulldata[fulldata$out==0,]
  tmpdata = fulldata
  

  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==1 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_1, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==1 & tmpdata$A_0==tmpdata$L1_0,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3_1 = (meany3tmp)

  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==1 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_1, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==1 & tmpdata$A_0==0,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3_2 = (meany3tmp)
  
  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==1 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_1 + L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==1 & tmpdata$A_0==tmpdata$L1_0,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3_3 = (meany3tmp)
  
  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==1 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_1, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==1 & tmpdata$A_0==tmpdata$L1_0,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3_4 = (meany3tmp)
  
  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==1 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_1 + L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==1 & tmpdata$A_0==0,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3_5 = (meany3tmp)
  
  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==1 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_1, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==1 & tmpdata$A_0==0,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3_6 = (meany3tmp)
  
  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==1 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_1 + L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==1 & tmpdata$A_0==0,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3_7 = (meany3tmp)
  
  ##################
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 & tmpdata$Y_2==1 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y3dat$y3pred = predict(y3fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_1 + L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$Y_1==1 & tmpdata$A_0==tmpdata$L1_0,]; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y2fito, newdata = y3dat, type="response"); 
  
  y3fit = glm(y3pred ~ L1_0, family = binomial(), data = y3dat) ; 
  y3dat = tmpdata; 
  y3dat$y3pred = predict(y3fit, newdata = y3dat, type="response")*predict(y1fito, newdata = y3dat, type="response"); 
  
  meany3tmp = c(mean(y3dat$y3pred))
  
  meany3_8 = (meany3tmp)
  
  
  

w1 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
  nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
  nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) 

w2 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
  nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
  nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) 

w3 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
  nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
  nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) 

w4 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
  nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
  nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) 

w5 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
  nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
  nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) 

w6 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
  nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
  nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) 

w7 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
  nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
  nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,])

w8 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
  nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
  nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,])


weights = c(w1,w2,w3,w4,w5,w6,w7,w8)
wsum = sum(weights)
myparam = 1-sum(c(meany3_1,meany3_2,meany3_3,meany3_4,meany3_5,meany3_6,meany3_7,meany3_8)*weights)
#proc.time() - ptm

return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"seq_reg_withextraeq_strat_3.csv")

stopCluster(cl)