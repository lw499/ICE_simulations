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
  library(geepack);library(MASS);
  library(dplyr)
  library(data.table)
  
  #library(devtools)
  #suppressPackageStartupMessages(install("gformula"))
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
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_3 &  
                    tmpdata$A_4==tmpdata$L1_4,]; 
  y5fito = glm(Y_5 ~ L1_4, family = binomial(), data = y5dat) ;

  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2 &  tmpdata$A_3==tmpdata$L1_3,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  L0fito = glm(L1_0 ~ 1, family = binomial(), data = tmpdata)
  
  ##monte carlo simulation
  #alpha0=summary(treatmodel)$coef[1,1]; alpha1=summary(treatmodel)$coef[2,1]; alpha2=0; alpha3=0;
  #beta0=summary(covmodel)$coef[1,1]; beta1=summary(covmodel)$coef[2,1]; beta2=0; beta3=0;
  #theta0=summary(ymodel)$coef[1,1]; theta1=summary(ymodel)$coef[2,1]; theta2=summary(ymodel)$coef[3,1]; theta3=0; theta4=0;

  #var = sqrt(mean((covmodel$y - fitted(covmodel))^2))
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4 <- L14 <- A4 <- Y5 <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)

  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  Y2 = plogis((summary(y2fito)$coef[1,1]+summary(y2fito)$coef[2,1]*L11))

  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  Y3 = plogis((summary(y3fito)$coef[1,1]+summary(y3fito)$coef[2,1]*L12))

  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  Y4 = plogis((summary(y4fito)$coef[1,1]+summary(y4fito)$coef[2,1]*L13))

  L14 = as.numeric(rexpit(rep(summary(L4fito)$coef[1,1],N)) | L13)
  Y5 = plogis((summary(y5fito)$coef[1,1]+summary(y5fito)$coef[2,1]*L14))
  
  #tmpdata = data.frame(L10, A0, Y1, L11, A1, Y2, L12, A2, Y3, L13, A3, Y4, L14, A4, Y5)
  #tmpdata = data.frame(Y1, Y2, Y3, Y4, Y5)
  
  risk1 = mean(Y1)
  risk2 = 1-mean((1-Y1)*(1-Y2))
  risk3 = 1-mean((1-Y1)*(1-Y2)*(1-Y3))
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  risk5 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4)*(1-Y5))
  
  myparam = c(risk1, risk2, risk3, risk4, risk5); #myparam
  #proc.time() - ptm
  
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"gform_sim_strat.csv")

stopCluster(cl)
