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
  #library(gformula)
  #library(Hmisc)
  
  #library(reshape2)  #do not use for data frame only
  rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
  
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
  dffull[, paste("lag_L") := shift(L1, 1, NA, type='lag'), by=id]
  #dffull[, paste("lag.L1",i, sep="") :=shift(L1, i, NA,type='lag'),by=id]
  #dffull[, paste("lag.L2",i, sep="") :=shift(L2, i, NA, type='lag'),by=id]
  
  covmodel<- glm(L1 ~ lag_A, data = dffull[!is.na(dffull$lag_L) & dffull$lag_L==0,], family=binomial())
  treatmodel <- glm(A ~ L1, data = dffull[!is.na(dffull$lag_A) & dffull$lag_A==0,], family=binomial())
  ymodel <- glm(Y ~ A + L1, data = dffull, family=binomial())

  ##monte carlo simulation
  alpha0=summary(treatmodel)$coef[1,1]; alpha1=summary(treatmodel)$coef[2,1]; alpha2=0; alpha3=0;
  beta0=summary(covmodel)$coef[1,1]; beta1=summary(covmodel)$coef[2,1]; beta2=0; beta3=0;
  theta0=summary(ymodel)$coef[1,1]; theta1=summary(ymodel)$coef[2,1]; theta2=summary(ymodel)$coef[3,1]; theta3=0; theta4=0;

  m = 1 #grace period is 2, plus one for Lk measured either at the end of month or beginning of month
  N=10000
  ua <- rep(TRUE, N)
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4 <- L14 <- A4 <- Y5 <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  A0 <- rexpit(alpha0+alpha1*L10)
  A0 = ifelse(L10==0, 0, A0)
  mtmp = ifelse(L10==1, 0, NA) # this keeps track of when threshold is crossed
  Y1 <- plogis(theta0+theta1*A0+theta2*L10)
  ua <- ua
  
  L11[ua] = as.numeric(rexpit((beta0+beta1*A0+beta2*L10)[ua]) | L10[ua])
  A1[ua] = as.numeric(rexpit((alpha0+alpha1*L11)[ua]) | A0[ua])
  A1[ua] = ifelse(L11[ua]==0, 0, A1[ua]) # If L11 is 0 (threshold is not crossed), then A1 is 0
  mtmp[ua] = ifelse(is.na(mtmp[ua]) & L11[ua]==1, 1, mtmp[ua])
  A1[ua] = ifelse(!is.na(mtmp[ua]) & mtmp[ua]+m==1, 1, A1[ua])
  Y2[ua] = plogis((theta0+theta1*A1+theta2*L11)[ua])
  ua <- ua
  
  L12[ua] = as.numeric(rexpit((beta0+beta1*A1+beta2*L11)[ua]) | L11[ua])
  A2[ua] = as.numeric(rexpit((alpha0+alpha1*L12)[ua]) | A1[ua])
  A2[ua] = ifelse(L12[ua]==0,0,A2[ua])
  mtmp[ua] = ifelse(is.na(mtmp[ua]) & L12[ua]==1, 2, mtmp[ua])
  A2[ua] = ifelse(!is.na(mtmp[ua]) & mtmp[ua]+m==2, 1, A2[ua])
  Y3[ua] = plogis((theta0+theta1*A2+theta2*L12)[ua])
  ua <- ua
  
  L13[ua] = as.numeric(rexpit((beta0+beta1*A2+beta2*L12)[ua]) | L12[ua])
  A3[ua] = as.numeric(rexpit((alpha0+alpha1*L13)[ua]) | A2[ua])
  A3[ua] = ifelse(L13[ua]==0,0,A3[ua])
  mtmp[ua] = ifelse(is.na(mtmp[ua]) & L13[ua]==1, 3, mtmp[ua])
  A3[ua] = ifelse(!is.na(mtmp[ua]) & mtmp[ua]+m==3, 1, A3[ua])
  Y4[ua] = plogis((theta0+theta1*A3+theta2*L13)[ua])
  ua <- ua
  
  L14[ua] = as.numeric(rexpit((beta0+beta1*A3+beta2*L13)[ua]) | L13[ua])
  A4[ua] = as.numeric(rexpit((alpha0+alpha1*L14)[ua]) | A3[ua])
  A4[ua] = ifelse(L14[ua]==0,0,A4[ua])
  mtmp[ua] = ifelse(is.na(mtmp[ua]) & L14[ua]==1, 4, mtmp[ua])
  A4[ua] = ifelse(!is.na(mtmp[ua]) & mtmp[ua]+m==4, 1, A4[ua])
  Y5[ua] = plogis((theta0+theta1*A4+theta2*L14)[ua])
  
  #tmpdata = data.frame(L10, A0, Y1, L11, A1, Y2, L12, A2, Y3, L13, A3, Y4, L14, A4, Y5)
  risk1 = mean(Y1)
  risk2 = 1-mean((1-Y1)*(1-Y2))
  risk3 = 1-mean((1-Y1)*(1-Y2)*(1-Y3))
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  risk5 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4)*(1-Y5))
  
  myparam = c(risk1, risk2, risk3, risk4, risk5)
  
    return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"gform_sim.csv")

stopCluster(cl)
