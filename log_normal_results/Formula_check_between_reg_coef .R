rho = 0.6
Sigma = matrix(c(1,rho,rho,1), 2, 2)
mu = c(0,0)
no_of_sim = 1000
beta_1 = matrix(NA, no_of_sim, 3)
for(i in 1:no_of_sim)
{
  X = mvrnorm(1000, mu, Sigma)
  Y = X%*%c(log(10.3),log(4.3)) + rnorm(1000,0,1)
  theta_1 = as.numeric(lm(Y ~ -1 + X[,1])$coeff)
  theta_2 = as.numeric(lm(Y ~ -1 + X[,2])$coeff)
  beta_1[i,1] = (theta_1 - rho*theta_2)/(1-rho^2)
  beta_1[i,2] = as.numeric(lm(Y ~ -1 + X[,1] + X[,2])$coeff[1])
  beta_1[i,3] = theta_1 - rho* as.numeric(lm(Y ~ -1 + X[,1] + X[,2])$coeff[2])
}

