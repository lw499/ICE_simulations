#library(dplyr)
#library(data.table)
#library(parallel)

#startTime <-proc.time()
#cl <- makeCluster(detectCores()-1)
set.seed(5678)

rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
N <- 10000000
Ndat <- 5
alpha0=-1.5; alpha1=1; alpha2=0; alpha3=0;
beta0=-2; beta1=-2; beta2=0; beta3=0;
theta0=-2; theta1=-2; theta2=1; theta3=0; theta4=0;
eta0=-3; eta1=-1; eta2=0.75; eta3=0; eta4=0;
gamma0=-2; gamma1=-1; gamma2=-1; gamma3=0; gamma4=0
sigma=0.1

m = 1 #grace period

ua <- rep(TRUE, N)
L10 <- L20 <- A0 <- Y1 <- L11 <- L21 <- A1 <- Y2 <- L12 <- L22 <- A2 <- Y3 <- L13 <- L23 <- A3 <- Y4 <- L14 <- L24 <- A4 <- Y5 <- L15 <- L25 <- A5 <- Y6 <- L16 <- L26 <- A6 <- Y7 <- as.numeric(rep(NA, N))
ids <- as.list(1:N)
L10 <- rbinom(N, 1, 0.5)
L20 <- rnorm(N, mean=beta0, sd=sigma)
A0 <- rexpit(alpha0+alpha1*L10+alpha2*L20)
A0 = ifelse(L10==0, 0, A0)
    mtmp = ifelse(L10==1, 0, NA) # this keeps track of when threshold is crossed
Y1 <- rbinom(N, 1, plogis(theta0+theta1*A0+theta2*L10+theta3*L20))
ua <- ua & !Y1

L11[ua] = as.numeric(rexpit((beta0+beta1*A0+beta2*L10)[ua]) | L10[ua])
L21[ua] = rnorm(sum(ua), mean=(beta0+beta1*A0+beta2*L10+beta3*L20)[ua], sd=sigma)
A1[ua] = as.numeric(rexpit((alpha0+alpha1*L11+alpha2*L21+alpha3*1)[ua]) | A0[ua])
    A1[ua] = ifelse(L11[ua]==0, 0, A1[ua]) # If L11 is 0 (threshold is not crossed), then A1 is 0
    mtmp[ua] = ifelse(is.na(mtmp[ua]) & L11[ua]==1, 1, mtmp[ua])
    A1[ua] = ifelse(!is.na(mtmp[ua]) & mtmp[ua]+m==1, 1, A1[ua])
Y2[ua] = rexpit((theta0+theta1*A1+theta2*L11+theta3*L21+theta4*1)[ua])
ua <- ua & !Y2

L12[ua] = as.numeric(rexpit((beta0+beta1*A1+beta2*L11)[ua]) | L11[ua])
L22[ua] = rnorm(sum(ua), mean=(beta0+beta1*A1+beta2*L11+beta3*L21)[ua], sd=sigma)
A2[ua] = as.numeric(rexpit((alpha0+alpha1*L12+alpha2*L22+alpha3*2)[ua]) | A1[ua])
    A2[ua] = ifelse(L12[ua]==0,0,A2[ua])
    mtmp[ua] = ifelse(is.na(mtmp[ua]) & L12[ua]==1, 2, mtmp[ua])
    A2[ua] = ifelse(!is.na(mtmp[ua]) & mtmp[ua]+m==2, 1, A2[ua])
Y3[ua] = rexpit((theta0+theta1*A2+theta2*L12+theta3*L22+theta4*2)[ua])
ua <- ua & !Y3

L13[ua] = as.numeric(rexpit((beta0+beta1*A2+beta2*L12)[ua]) | L12[ua])
L23[ua] = rnorm(sum(ua), mean=(beta0+beta1*A2+beta2*L12+beta3*L22)[ua], sd=sigma)
A3[ua] = as.numeric(rexpit((alpha0+alpha1*L13+alpha2*L23+alpha3*3)[ua]) | A2[ua])
    A3[ua] = ifelse(L13[ua]==0,0,A3[ua])
    mtmp[ua] = ifelse(is.na(mtmp[ua]) & L13[ua]==1, 3, mtmp[ua])
    A3[ua] = ifelse(!is.na(mtmp[ua]) & mtmp[ua]+m==3, 1, A3[ua])
Y4[ua] = rexpit((theta0+theta1*A3+theta2*L13+theta3*L23+theta4*3)[ua])
ua <- ua & !Y4

L14[ua] = as.numeric(rexpit((beta0+beta1*A3+beta2*L13)[ua]) | L13[ua])
L24[ua] = rnorm(sum(ua), mean=(beta0+beta1*A3+beta2*L13+beta3*L23)[ua], sd=sigma)
A4[ua] = as.numeric(rexpit((alpha0+alpha1*L14+alpha2*L24+alpha3*4)[ua]) | A3[ua])
    A4[ua] = ifelse(L14[ua]==0,0,A4[ua])
    mtmp[ua] = ifelse(is.na(mtmp[ua]) & L14[ua]==1, 4, mtmp[ua])
    A4[ua] = ifelse(!is.na(mtmp[ua]) & mtmp[ua]+m==4, 1, A4[ua])
Y5[ua] = rexpit((theta0+theta1*A4+theta2*L14+theta3*L24+theta4*4)[ua])

tmpdata = data.frame(L10, L20, A0, Y1, L11, L21, A1, Y2, L12, L22, A2, Y3, L13, L23, A3, Y4, L14, L24, A4, Y5)
#tmpdata = data.frame(Y1, Y2, Y3, Y4, Y5, Y6, Y7)

tmpdata$Y2 = ifelse(tmpdata$Y1==1,1,tmpdata$Y2)
tmpdata$Y3 = ifelse(!is.na(tmpdata$Y2) & tmpdata$Y2==1,1,tmpdata$Y3)
tmpdata$Y4 = ifelse(!is.na(tmpdata$Y3) & tmpdata$Y3==1,1,tmpdata$Y4)
tmpdata$Y5 = ifelse(!is.na(tmpdata$Y4) & tmpdata$Y4==1,1,tmpdata$Y5)

c(mean(tmpdata$Y1),mean(tmpdata$Y2),mean(tmpdata$Y3),mean(tmpdata$Y4),mean(tmpdata$Y5))
#Random: 0.1523100 0.2275338 0.2914871 0.3461815 0.3936240
