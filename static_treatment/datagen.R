# Function to create survival data in long form
# It simulates ONE INDIVIDUAL's data at t = {0, 1, 2, ..., K-1} (outcomes realized at t = {1, 2, ..., K})
# Y is binary outcome
# Right censoring is present (C)
# A is binary treatment
# No competing risk
# 2 covariates (L1, L2)

# Inputs are i: id, df: data frame specified above, K: max number of time points
# Alphas, betas, thetas, etas, and sigma are user-provided parameters for data generating functions
datagen <- function(i, K, alpha0, alpha1, alpha2, alpha3, beta0, beta1, beta2, beta3,
                    theta0, theta1, theta2, theta3, theta4, eta0, eta1, eta2, eta3, eta4, sigma){
  id <- as.numeric(i)
  # Data at time 0
  L1 <- rbinom(1, 1, 0.5)
  L2 <- rnorm(1, mean=beta0, sd=sigma)
  A <- rbinom(1, 1, plogis(alpha0+alpha1*L1+alpha2*L2))
  C <- rbinom(1, 1, plogis(eta0+eta1*A+eta2*L1+eta3*L2))
  Y <- rbinom(1, 1, plogis(theta0+theta1*A+theta2*L1+theta3*L2))
  temp <- data.frame(id = id, t0 = 0, L1, L2, A, C, Y)
  
  if (temp$C==1){
    temp$Y <- NA
  }else if (temp$Y!=1){
    for (j in 2:K){
      t0 <- j-1
      L1star <- as.numeric(rbinom(1, 1, plogis(beta0+beta1*temp$A[j-1]+beta2*temp$L1[j-1])) | temp$L1[j-1])
      L2star <- rnorm(1, mean=beta0+beta1*temp$A[j-1]+beta2*temp$L1[j-1]+beta3*temp$L2[j-1], sd=sigma)
      temp_L1 <- c(temp$L1, L1star)
      #Astar <- rbinom(1, 1, plogis(alpha0+alpha1*cumsum(temp_L1)[j]/seq_along(temp_L1)[j]
      #                             +alpha2*L2star+alpha3*t0))
      Astar <- as.numeric(rbinom(1, 1, plogis(alpha0+alpha1*L1star+alpha2*L2star+alpha3*t0)) | temp$A[j-1])
      temp_A <- c(temp$A, Astar)
      #Cstar <- rbinom(1, 1, plogis(eta0+eta1*cumsum(temp_A)[j]/seq_along(temp_A)[j]
      #                             +eta2*cumsum(temp_L1)[j]/seq_along(temp_L1)[j]
      #                             +eta3*L2star+eta4*t0))
      Cstar <- rbinom(1, 1, plogis(eta0+eta1*Astar+eta2*L1star+eta3*L2star+eta4*t0))
      if (Cstar==1){
        Ystar <- NA
        temp <- rbind(temp, c(id, t0, L1star, L2star, Astar, Cstar, Ystar))
        break
      }
      else{
        #Ystar <- rbinom(1, 1, plogis(theta0+theta1*cumsum(temp_A)[j]/seq_along(temp_A)[j]
        #                             +theta2*cumsum(temp_L1)[j]/seq_along(temp_L1)[j]
        #                             +theta3*L2star+theta4*t0))
        Ystar <- rbinom(1, 1, plogis(theta0+theta1*Astar+theta2*L1star+theta3*L2star+theta4*t0))
      }
      temp <- rbind(temp, c(id, t0, L1star, L2star, Astar, Cstar, Ystar))
      if(Ystar==1){break}
    }
  }
  return (temp)
}
