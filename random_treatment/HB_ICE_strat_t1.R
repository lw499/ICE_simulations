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
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y1dat = tmpdata; 
  y1dat$y1pred = predict(y1fito, newdata = y1dat, type="response"); 
  
  meany1tmp = c(mean(y1dat$y1pred))
  
  meany1_1 = (meany1tmp)
  
  ##################
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  
  y1dat = tmpdata; 
  y1dat$y1pred = predict(y1fito, newdata = y1dat, type="response"); 
  
  meany1tmp = c(mean(y1dat$y1pred))
  
  meany1_2 = (meany1tmp)


w1 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,])
w2 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,])


weights = c(w1,w2)
wsum = sum(weights)
myparam = 1-sum(c(meany1_1,meany1_2)*weights)
#proc.time() - ptm

return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"seq_reg_withextraeq_strat_1.csv")

stopCluster(cl)