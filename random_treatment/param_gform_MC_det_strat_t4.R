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
  library(geepack);library(MASS);
  library(dplyr)
  library(data.table)
  setDTthreads(1)
  
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
  
  #####
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

  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")

  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")

  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))

  meany4_1 = risk4

  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2 &  tmpdata$A_3==tmpdata$L1_3,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_2 = risk4
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2 &  tmpdata$A_3==tmpdata$L1_2,]; 
  y4fito = glm(Y_4 ~ L1_3 + L1_2, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_3 = risk4 ###
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2 &  tmpdata$A_3==tmpdata$L1_2,]; 
  y4fito = glm(Y_4 ~ L1_3 + L1_2, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_4 = risk4

  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2 &  tmpdata$A_3==tmpdata$L1_3,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_5 = risk4
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1 &  tmpdata$A_3==tmpdata$L1_2,]; 
  y4fito = glm(Y_4 ~ L1_3 + L1_2, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_6 = risk4
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1 &  tmpdata$A_3==tmpdata$L1_3,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1,family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_0,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_7 = risk4 
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1 &  tmpdata$A_3==tmpdata$L1_2,]; 
  y4fito = glm(Y_4 ~ L1_3 + L1_2, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_8 = risk4 
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1 &  tmpdata$A_3==tmpdata$L1_3,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==tmpdata$L1_0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_9 = risk4 
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1 &  tmpdata$A_3==tmpdata$L1_2,]; 
  y4fito = glm(Y_4 ~ L1_3 + L1_2, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 ,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_10 = risk4 
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1 &  tmpdata$A_3==tmpdata$L1_3,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 ,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_11 = risk4 
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1 &  tmpdata$A_3==tmpdata$L1_2,]; 
  y4fito = glm(Y_4 ~ L1_3 + L1_2, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 ,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_12 = risk4 
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1 &  tmpdata$A_3==tmpdata$L1_3,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_1,]; 
  y3fito = glm(Y_3 ~ L1_2 + L1_1, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_1 ,]; 
  y2fito = glm(Y_2 ~ L1_1, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_13 = risk4 
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2 &  tmpdata$A_3==tmpdata$L1_2,]; 
  y4fito = glm(Y_4 ~ L1_3 + L1_2, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 ,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_14 = risk4 
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2 &  tmpdata$A_3==tmpdata$L1_3,]; 
  y4fito = glm(Y_4 ~ L1_3, family = binomial(), data = y4dat) ; 
  L4fito = glm(L1_4 ~ 1, family = binomial(), data = y4dat[y4dat$L1_3==0,]) ; 
  
  y3dat = tmpdata[!is.na(tmpdata$C_3) & tmpdata$C_3==0 & 
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 &  
                    tmpdata$A_2==tmpdata$L1_2,]; 
  y3fito = glm(Y_3 ~ L1_2, family = binomial(), data = y3dat) ; 
  L3fito = glm(L1_3 ~ 1, family = binomial(), data = y3dat[y3dat$L1_2==0,]) ; 
  
  y2dat = tmpdata[!is.na(tmpdata$C_2) & tmpdata$C_2==0 &
                    tmpdata$A_0==0 & tmpdata$A_1==tmpdata$L1_0 ,]; 
  y2fito = glm(Y_2 ~ L1_1 + L1_0, family = binomial(), data = y2dat) ; 
  L2fito = glm(L1_2 ~ 1, family = binomial(), data = y2dat[y2dat$L1_1==0,]) ; 
  
  y1dat = tmpdata[!is.na(tmpdata$C_1) & tmpdata$C_1==0 & tmpdata$A_0==0,]; 
  y1fito = glm(Y_1 ~ L1_0, family = binomial(), data = y1dat) ; 
  L1fito = glm(L1_1 ~ 1, family = binomial(), data = y1dat[y1dat$L1_0==0,]) ; 
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_15 = risk4 
  
  
  #####
  y4dat = tmpdata[!is.na(tmpdata$C_4) & tmpdata$C_4==0 &
                    tmpdata$A_0==tmpdata$L1_0 & tmpdata$A_1==tmpdata$L1_1 &  
                    tmpdata$A_2==tmpdata$L1_2 &  tmpdata$A_3==tmpdata$L1_2,]; 
  y4fito = glm(Y_4 ~ L1_3 + L1_2, family = binomial(), data = y4dat) ; 
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
  
  N=10000
  L10 <- A0 <- Y1 <- L11 <- A1 <- Y2 <- L12 <- A2 <- Y3 <- L13 <- A3 <- Y4  <- as.numeric(rep(NA, N))
  ids <- as.list(1:N)
  L10 <- sample(dffull[dffull$t0==0,]$L1, N, replace=T)
  Y1 <- plogis(summary(y1fito)$coef[1,1]+summary(y1fito)$coef[2,1]*L10)
  
  L11 = as.numeric(rexpit(rep(summary(L1fito)$coef[1,1],N)) | L10)
  L1_1 = L11; L1_0 = L10; tmp = as.data.frame(cbind(L1_1,L1_0))
  Y2 = predict(y2fito, newdata=tmp, type="response")
  
  L12 = as.numeric(rexpit(rep(summary(L2fito)$coef[1,1],N)) | L11)
  L1_2 = L12; L1_1 = L11; tmp = as.data.frame(cbind(L1_2,L1_1))
  Y3 = predict(y3fito, newdata=tmp, type="response")
  
  L13 = as.numeric(rexpit(rep(summary(L3fito)$coef[1,1],N)) | L12)
  L1_3 = L13; L1_2 = L12; tmp = as.data.frame(cbind(L1_3,L1_2))
  Y4 = predict(y4fito, newdata=tmp, type="response")
  risk4 = 1-mean((1-Y1)*(1-Y2)*(1-Y3)*(1-Y4))
  
  meany4_16 = risk4 
  
  
  w1 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_3,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])
  
  w2 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_3,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])
  
  
  w3 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,]) 
  
  w4 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])
  
  w5 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_3,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,]) 
  
  w6 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1 & tmpdata$A_3==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,]) 
  
  w7 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1 & tmpdata$A_3==tmpdata$L1_3,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,]) 
  
  w8 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1 & tmpdata$A_3==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])
  
  w9 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1 & tmpdata$A_3==tmpdata$L1_3,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])
  
  w10 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1 & tmpdata$A_3==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,])
  
  w11 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1 & tmpdata$A_3==tmpdata$L1_3,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_1,])
  
  w12 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1 & tmpdata$A_3==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])
  
  w13 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1 & tmpdata$A_3==tmpdata$L1_3,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_1,])
  
  w14 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])
  
  w15 = nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_0,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_0 ,]) *
    nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_3,])/nrow(tmpdata[tmpdata$A_0==0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_0 & tmpdata$A_2==tmpdata$L1_2,])
  
  w16 = nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==1,])/nrow(tmpdata[tmpdata$L1_0==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1 & tmpdata$A_1==tmpdata$L1_1,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==1,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==1 & tmpdata$A_1==tmpdata$L1_1 ,]) *
    nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2 & tmpdata$A_3==tmpdata$L1_2,])/nrow(tmpdata[tmpdata$A_0==tmpdata$L1_0 & tmpdata$L1_0==0 & tmpdata$L1_1==0 & tmpdata$L1_2==0 & tmpdata$L1_3==1 & tmpdata$A_1==tmpdata$L1_1 & tmpdata$A_2==tmpdata$L1_2,])
  
  
  
  
  
  weights = c(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16)
  wsum = sum(weights)
  myparam = sum(c(meany4_1,meany4_2,meany4_3,meany4_4,meany4_5,meany4_6,meany4_7,meany4_8,meany4_9,meany4_10,meany4_11,meany4_12,meany4_13,meany4_14,meany4_15,meany4_16)*weights)
  # [1] 0.7353390 0.6322651 0.6236331 0.7126780 0.7198491 0.7045067 0.7117274 0.7199476 0.7271625 0.6004848 0.6081815 0.6174468 0.6250799 0.6076647 0.6152843 0.7281718
  
  
  
  
  
  
    
  return(myparam)
}
test = foreach(m=1:1000) %dopar% myfunc(m)
test2 = do.call("rbind", test)

write.csv(test2,"gform_sim_strat_4.csv")

stopCluster(cl)
