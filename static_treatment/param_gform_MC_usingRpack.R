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
  
  #library(devtools)
  #suppressPackageStartupMessages(install("gformula"))
  library(gformula)
  #library(Hmisc)
  
  #library(reshape2)  #do not use for data frame only
  
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
  #dffull[, paste("lag.Y",i, sep="") := shift(Y, i, NA, type='lag'), by=id]
  #dffull[, paste("lag.C",i, sep="") := shift(C, i, NA, type='lag'), by=id]
  dffull[, paste("lag_A") := shift(A, 1, NA, type='lag'), by=id]
  dffull[, paste("lag_L1") := shift(L1, 1, NA, type='lag'), by=id]
  #dffull[, paste("lag.L1",i, sep="") :=shift(L1, i, NA,type='lag'),by=id]
  #dffull[, paste("lag.L2",i, sep="") :=shift(L2, i, NA, type='lag'),by=id]
  
  #ptm <- proc.time()
  
  time_points <- 5
  covnames <- c('L1', 'L2', 'A')
  covtypes <- c('binary', 'bounded normal', 'binary')
  covparams <- list(covlink = c('logit', 'identity', 'logit'),
                    covmodels = c(L1 ~ lag_A,
                                  L2 ~ lag_A,
                                  A ~ L1))
  ymodel <- Y ~ L1 + A 
  intvars <- c('A', 'A')
  interventions <- list(c(static, 1),c(static, 0))
  int_descript <- c('Always treat', 'Never treat')
  histories <- c(lagged)
  histvars <- c('L1', 'L2', 'A')
  #nsimul = 10000
  restrictions = list(c('L1', 'lag_L1 == 0', simple_restriction, 1))
  suppressMessages(gform_basic <- gformula_survival(obs_data = dffull, time_points = time_points,
                                               covnames = covnames, covtypes = covtypes,
                                               covparams = covparams, ymodel = ymodel,
                                               intvars = intvars, interventions = interventions,
                                               restrictions = restrictions,nsimul = 10000,
                                               int_descript = int_descript,
                                               histories = histories, histvars = histvars,
                                               seed = 1234))
  myparam = rbind(c(1, gform_basic$`Result table`[Interv.==1,]$`Est. risk`), c(2, gform_basic$`Result table`[Interv.==2,]$`Est. risk`))
  #proc.time() - ptm
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"gform_sim.csv")

stopCluster(cl)
