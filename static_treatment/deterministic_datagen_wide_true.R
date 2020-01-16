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

m = 2 #grace period is 2, plus one for Lk measured either at the end of month or beginning of month

ua <- rep(TRUE, N)
L10 <- L20 <- A0 <- Y1 <- L11 <- L21 <- A1 <- Y2 <- L12 <- L22 <- A2 <- Y3 <- L13 <- L23 <- A3 <- Y4 <- L14 <- L24 <- A4 <- Y5 <- L15 <- L25 <- A5 <- Y6 <- L16 <- L26 <- A6 <- Y7 <- as.numeric(rep(NA, N))
ids <- as.list(1:N)
L10 <- rbinom(N, 1, 0.5)
L20 <- rnorm(N, mean=beta0, sd=sigma)
A0 <- rep(1, N)
Y1 <- rbinom(N, 1, plogis(theta0+theta1*A0+theta2*L10+theta3*L20))
ua <- ua & !Y1

L11[ua] = as.numeric(rexpit((beta0+beta1*A0+beta2*L10)[ua]) | L10[ua])
L21[ua] = rnorm(sum(ua), mean=(beta0+beta1*A0+beta2*L10+beta3*L20)[ua], sd=sigma)
A1[ua] = rep(1,N)[ua]
#A1[ua] = L11[ua] # If L11 is 0 (threshold is not crossed), then A1 is 0
Y2[ua] = rexpit((theta0+theta1*A1+theta2*L11+theta3*L21+theta4*1)[ua])
ua <- ua & !Y2

L12[ua] = as.numeric(rexpit((beta0+beta1*A1+beta2*L11)[ua]) | L11[ua])
L22[ua] = rnorm(sum(ua), mean=(beta0+beta1*A1+beta2*L11+beta3*L21)[ua], sd=sigma)
A2[ua] = rep(1, N)[ua]
#A2[ua] = L12[ua]
Y3[ua] = rexpit((theta0+theta1*A2+theta2*L12+theta3*L22+theta4*2)[ua])
ua <- ua & !Y3

L13[ua] = as.numeric(rexpit((beta0+beta1*A2+beta2*L12)[ua]) | L12[ua])
L23[ua] = rnorm(sum(ua), mean=(beta0+beta1*A2+beta2*L12+beta3*L22)[ua], sd=sigma)
A3[ua] = rep(1, N)[ua]
#A3[ua] = L13[ua]
Y4[ua] = rexpit((theta0+theta1*A3+theta2*L13+theta3*L23+theta4*3)[ua])
ua <- ua & !Y4

L14[ua] = as.numeric(rexpit((beta0+beta1*A3+beta2*L13)[ua]) | L13[ua])
L24[ua] = rnorm(sum(ua), mean=(beta0+beta1*A3+beta2*L13+beta3*L23)[ua], sd=sigma)
A4[ua] = rep(1, N)[ua]
#A4[ua] = L14[ua]
Y5[ua] = rexpit((theta0+theta1*A4+theta2*L14+theta3*L24+theta4*4)[ua])

tmpdata = data.frame(L10, L20, A0, Y1, L11, L21, A1, Y2, L12, L22, A2, Y3, L13, L23, A3, Y4, L14, L24, A4, Y5)
#tmpdata = data.frame(Y1, Y2, Y3, Y4, Y5, Y6, Y7)

tmpdata$Y2 = ifelse(tmpdata$Y1==1,1,tmpdata$Y2)
tmpdata$Y3 = ifelse(!is.na(tmpdata$Y2) & tmpdata$Y2==1,1,tmpdata$Y3)
tmpdata$Y4 = ifelse(!is.na(tmpdata$Y3) & tmpdata$Y3==1,1,tmpdata$Y4)
tmpdata$Y5 = ifelse(!is.na(tmpdata$Y4) & tmpdata$Y4==1,1,tmpdata$Y5)


c(mean(tmpdata$Y1),mean(tmpdata$Y2),mean(tmpdata$Y3),mean(tmpdata$Y4),mean(tmpdata$Y5))
#deterministic:  0.0326666 0.0643971 0.0950325 0.1247307 0.1535464
#deterministic: 0.1940054 0.3527089 0.4817619 0.5861110 0.6704213

#(nrow(datatrue[datatrue$Y2==1 & datatrue$Y1==0,])/nrow(datatrue[datatrue$Y1==0,]))
#0.1934693
