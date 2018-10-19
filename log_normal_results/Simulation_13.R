#Inclusion Criterion on reference dataset

library(MASS)
library(stats)
#####
### Basic setting
#####
d.X = 3 # number of covariates.
mu.1 = matrix(rep(0,d.X), nrow=d.X) # mean vector of the covariates.
mu.2 = matrix(rep(0,d.X), nrow=d.X) # mean vector of the covariates.
mu.3 = matrix(rep(0,d.X), nrow=d.X) # mean vector of the covariates.

r1.1 = 0.3 # correlation coefficient of the covariates.
r2.1 = 0.6
r3.1 = 0.1
Sigma.1 = matrix( 
  c(1, r1.1, r2.1,  
    r1.1, 1, r3.1, 
    r2.1, r3.1, 1), 
  nrow=d.X, 
  ncol=d.X) # covariance matrix of the covariates.

r1.2 = 0.6 # correlation coefficient of the covariates.
r2.2 = 1.2
r3.2 = 0.2
Sigma.2 = matrix( 
  c(2, r1.2, r2.2,  
    r1.2, 2, r3.2, 
    r2.2, r3.2, 2), 
  nrow=d.X, 
  ncol=d.X) # covariance matrix of the covariates.


r1.3 = 0.15 # correlation coefficient of the covariates.
r2.3 = 0.3
r3.3 = 0.05
Sigma.3 = matrix( 
  c(0.5, r1.3, r2.3,  
    r1.3, 0.5, r3.3, 
    r2.3, r3.3, 0.5), 
  nrow=d.X, 
  ncol=d.X) # covariance matrix of the covariates.


beta.star = matrix(c(-1.2, log(1.3), log(1.3), log(1.3)),nrow = d.X+1) # beta.star
#beta.star = matrix(c(-1.2, 0.26, 0.26, 0.26),nrow = d.X+1) # beta.star
#beta.star = matrix(c(-3, 1, 2, 3),nrow = d.X+1) # beta.star


corr.data.m1 = list()
corr.data.m2 = list()
corr.data.m3 = list()
corr.data.X.rf = list()


n1 = 600 # sample size of the 1st data set.
n2 = 500 # 2nd
n3 = 2000 # 3rd

n = 100
no.of.simulations = 1000
sim.matrix = matrix(NA, 1000, 8)
asym_var_simu = list()
sim = 0
for(sim in 1:no.of.simulations)
{
  set.seed(sim)
  # Generate the reference data set
  X.rf = mvrnorm(n = n, mu.1, Sigma.1)
  X.rf = X.rf[X.rf[,1]>-0.5 & X.rf[,2]<0.5, ]
  
  
  # Generate data set 1. m1 means model 1.
  X.m1 = mvrnorm(n = n1, mu.1, Sigma.1) # Generate the covariates.
  #X.m1 = X.m1[X.m1[,1] > -0.5 & X.m1[,2] < 0.5, ]
  X.m1.1 = cbind(rep(1, dim(X.m1)[1]), X.m1) # Add a column of 1's to X.m1.
  p.m1 = 1/(1+exp(-X.m1.1%*%beta.star)) # the vector of probabilities
  Y.m1 = rbinom(length(p.m1), size=1, p.m1) # the Bernoulli responses
  # print(p.m1[1])
  # print(mean(Y.m1))
  # print(mean(p.m1))
  
  # Generate data set 2. m1 means model 2.
  X.m2 = mvrnorm(n = n2, mu.2, Sigma.1)
  X.m2.1 = cbind(rep(1, n2), X.m2)
  p.m2 = 1/(1+exp(-X.m2.1%*%beta.star))
  Y.m2 = rbinom(n2, size=1, p.m2)
  
  # Generate data set 3. m1 means model 3.
  X.m3 = mvrnorm(n = n3, mu.3, Sigma.1)
  #X.m3 = X.m3[X.m3[,1]>0, ]
  X.m3.1 = cbind(rep(1, dim(X.m3)[1]), X.m3)
  p.m3 = 1/(1+exp(-X.m3.1%*%beta.star))
  Y.m3 = rbinom(length(p.m3), size=1, p.m3)
  
  #####
  ### Create data sets in the format of data frame.
  #####
  data.m1 = data.frame(Y=Y.m1, X.m1)
  data.m1 = data.m1[data.m1$X1>-0.5 & data.m1$X2<0.5, ]
  data.m2 = data.frame(Y=Y.m2, X.m2)
  data.m3 = data.frame(Y=Y.m3, X.m3)
  data.m3 = data.m3[data.m3$X1>0, ]
  # str(data.m1)
  
  
  #---Sample correlation matrix---#
  corr.data.m1[[sim]] = cor(data.m1[,-1])
  corr.data.m2[[sim]] = cor(data.m2[,-1])
  corr.data.m3[[sim]] = cor(data.m3[,-1])
  corr.data.X.rf[[sim]] = cor(X.rf)
  
  #####
  ### Apply logistic regression with reduced models to the data sets 
  #####
  logit.m1 <- glm(Y ~ X1 + X2, data = data.m1, family = "binomial")
  # print(logit.m1)
  if(logit.m1$converged == FALSE)
  {
    print("glm for logit.m1 is not convergent.")
    next
  }
  
  logit.m2 <- glm(Y ~ X2 + X3, data = data.m2, family = "binomial")
  # print(logit.m2)
  if(logit.m2$converged == FALSE)
  {
    print("glm for logit.m2 is not convergent.")
    next
  }
  
  logit.m3 <- glm(Y ~ X1 + X3, data = data.m3, family = "binomial")
  # print(logit.m3)
  if(logit.m3$converged == FALSE)
  {
    print("glm for logit.m3 is not convergent.")
    next
  }
  
  
  #####
  ### Obtain the estimators of the parameters in the reduced models. 
  #####
  theta.m1 = logit.m1$coefficients
  theta.m2 = logit.m2$coefficients
  theta.m3 = logit.m3$coefficients
  
  
  
  #####
  ### Find the covariance matrix estimators for the reduced models
  #####
  
  
  
  #####
  # Basic notations for inputs
  #####
  
  K = 3 # Number of data sets
  
  A1 = c(1, 2) # index set A1, the indexes of the covariates of data set 1.
  A2 = c(2, 3) # index set A2
  A3 = c(1, 3) # index set A3
  
  
  X.m1.used = as.matrix(cbind(rep(1, dim(data.m1)[1]), data.m1[, A1+1, drop=F]))
  X.m2.used = as.matrix(cbind(rep(1, dim(data.m2)[1]), data.m2[, A2+1, drop=F]))
  X.m3.used = as.matrix(cbind(rep(1, dim(data.m3)[1]), data.m3[, A3+1, drop=F]))
  # str(X.m1.used)
  # str(X.m2.used)
  # str(X.m3.used)
  
  ##### Find Sigma.m1
  
  T.1 = matrix(rep(0, (length(A1)+1)^2), nrow=length(A1)+1)
  T.2 = T.1
  
  for (i in 1:dim(data.m1)[1])
  {
    a = as.vector(exp(-X.m1.used[i, , drop=F]%*%theta.m1))
    T.1 = T.1 + (a/(1+a)^2) * (t(X.m1.used[i, , drop=F])%*%X.m1.used[i, , drop=F])
  }
  
  for (i in 1:dim(data.m1)[1])
  {
    a = as.vector(1/( 1 + exp(-X.m1.used[i, , drop=F]%*%theta.m1)))
    T.2 = T.2 + (Y.m1[i]-a)^2 * (t(X.m1.used[i, , drop=F])%*%X.m1.used[i, , drop=F])
  }
  
  Sigma.m1 = solve(T.1)%*%T.2%*%solve(T.1) # This is actually Sigma.m1.n1. 
  
  ##### Find Sigma.m2
  
  T.1 = matrix(rep(0, (length(A2)+1)^2), nrow=length(A2)+1)
  T.2 = T.1
  
  for (i in 1:dim(data.m2)[1])
  {
    a = as.vector(exp(-X.m2.used[i, , drop=F]%*%theta.m2))
    T.1 = T.1 + (a/(1+a)^2) * (t(X.m2.used[i, , drop=F])%*%X.m2.used[i, , drop=F])
  }
  
  for (i in 1:dim(data.m2)[1])
  {
    a = as.vector(1/( 1 + exp(-X.m2.used[i, , drop=F]%*%theta.m2)))
    T.2 = T.2 + (Y.m2[i]-a)^2 * (t(X.m2.used[i, , drop=F])%*%X.m2.used[i, , drop=F])
  }
  
  Sigma.m2 = solve(T.1)%*%T.2%*%solve(T.1)
  
  
  ##### Find Sigma.m3
  
  T.1 = matrix(rep(0, (length(A3)+1)^2), nrow=length(A3)+1)
  T.2 = T.1
  
  for (i in 1:dim(data.m3)[1])
  {
    a = as.vector(exp(-X.m3.used[i, , drop=F]%*%theta.m3))
    T.1 = T.1 + (a/(1+a)^2) * (t(X.m3.used[i, , drop=F])%*%X.m3.used[i, , drop=F])
  }
  
  for (i in 1:dim(data.m3)[1])
  {
    a = as.vector(1/( 1 + exp(-X.m3.used[i, , drop=F]%*%theta.m3)))
    T.2 = T.2 + (Y.m3[i]-a)^2 * (t(X.m3.used[i, , drop=F])%*%X.m3.used[i, , drop=F])
  }
  
  Sigma.m3 = solve(T.1)%*%T.2%*%solve(T.1)
  
  
  
  
  
  # # Generate data set 1. m1 means model 1.
  # X.m1 = rmvnorm(n = n1, mu, Sigma) # Generate the covariates.
  # X.m1.1 = cbind(rep(1, n1), X.m1) # Add a column of 1's to X.m1.
  # p.m1 = 1/(1+exp(-X.m1.1%*%beta.star)) # the vector of probabilities
  # Y.m1 = rbinom(n1, size=1, p.m1) # the Bernoulli responses
  # # print(p.m1[1])
  # # print(mean(Y.m1))
  # # print(mean(p.m1))
  # 
  # # Generate data set 2. m1 means model 2.
  # X.m2 = rmvnorm(n = n2, mu, Sigma)
  # X.m2.1 = cbind(rep(1, n2), X.m2)
  # p.m2 = 1/(1+exp(-X.m2.1%*%beta.star))
  # Y.m2 = rbinom(n2, size=1, p.m2)
  # 
  # # Generate data set 3. m1 means model 3.
  # X.m3 = rmvnorm(n = n3, mu, Sigma)
  # X.m3.1 = cbind(rep(1, n3), X.m3)
  # p.m3 = 1/(1+exp(-X.m3.1%*%beta.star))
  # Y.m3 = rbinom(n3, size=1, p.m3)
  
  #####
  ### Create data sets in the format of data frame.
  #####
  # data.m1 = data.frame(Y=Y.full[1:n1], X.full[1:n1, ])
  # data.m2 = data.frame(Y=Y.full[((n1+1):(n1+n2))], X.full[((n1+1):(n1+n2)), ])
  # data.m3 = data.frame(Y=Y.full[((n1+n2+1):(n1+n2+n3))], X.full[((n1+n2+1):(n1+n2+n3)), ])
  # 
  # fit.1 = glm(Y ~ X2 + X3, family = binomial(), data = data.m1)
  # fit.2 = glm(Y ~ X3 + X4, family = binomial(), data = data.m2)
  # fit.3 = glm(Y ~ X2 + X4, family = binomial(), data = data.m3)
  
  names(theta.m1)=c("(Intercept)","Age","Height")
  names(theta.m2)=c("(Intercept)","Height", "Weight")
  names(theta.m3)=c("(Intercept)","Age", "Weight")
  
  ###now put in the GENMETA example
  
  study1 = list(Coeff=theta.m1,Covariance=Sigma.m1,Sample_size=n1)
  study2 = list(Coeff=theta.m2,Covariance=Sigma.m2,Sample_size=n2)
  study3 = list(Coeff=theta.m3,Covariance=Sigma.m3,Sample_size=n3)
  
  studies = list(study1,study2,study3)
  model = "logistic"
  
  reference = cbind(rep(1,dim(X.rf)[1]), X.rf)
  colnames(reference) = c("(Intercept)","Age","Height", "Weight")
  result.same = GENMETA(studies, reference, model, initial_val = c(-1.2, log(1.3), log(1.3), log(1.3)))
  
  asym_var_simu[[sim]] = result.same$Est.var.cov
  if(sum(is.na(result.same$Est.coeff)) == 0 && is.null(result.same$Est.var.cov) != TRUE)
  {
    sim.matrix[sim, ] = c(result.same$Est.coeff, diag(result.same$Est.var.cov))
  }
  
  print(sim)
}

Lower.CI = na.omit(sim.matrix)[,1:4] - 1.96*sqrt(na.omit(sim.matrix)[,5:8])
Upper.CI = na.omit(sim.matrix)[,1:4] + 1.96*sqrt(na.omit(sim.matrix)[,5:8])

Coverage.rate = matrix(0,dim(na.omit(sim.matrix))[1], 4)
for(k in 1:dim(na.omit(sim.matrix))[1])
{
  if(beta.star[1,1] >= Lower.CI[k,1] & beta.star[1,1] <= Upper.CI[k,1])
    Coverage.rate[k,1] = 1
  if(beta.star[2,1] >= Lower.CI[k,2] & beta.star[2,1] <= Upper.CI[k,2])
    Coverage.rate[k,2] = 1
  if(beta.star[3,1] >= Lower.CI[k,3] & beta.star[3,1] <= Upper.CI[k,3])
    Coverage.rate[k,3] = 1
  if(beta.star[4,1] >= Lower.CI[k,4] & beta.star[4,1] <= Upper.CI[k,4])
    Coverage.rate[k,4] = 1
}

Average.length = apply((Upper.CI-Lower.CI), 2, mean)
result = data.frame(cbind(apply(na.omit(sim.matrix), 2, mean)[1:4], beta.star, (apply(na.omit(sim.matrix), 2, mean)[1:4]- beta.star), sqrt(apply(na.omit(sim.matrix), 2, mean)[5:8]), sqrt(apply(na.omit(sim.matrix)[,1:4], 2, var))), apply(Coverage.rate,2, mean), Average.length)
colnames(result) = c("Coeff","True", "Bias", "ESE", "SE", "Coverage.rate", "AL")



write.csv(result, file = "Simulation_13.csv", row.names = F)

write.csv(Reduce("+", corr.data.m1) / length(corr.data.m1), file = "Simulation_13_study_1_corr.csv")
write.csv(Reduce("+", corr.data.m2) / length(corr.data.m2), file = "Simulation_13_study_2_corr.csv")
write.csv(Reduce("+", corr.data.m3) / length(corr.data.m3), file = "Simulation_13_study_3_corr.csv")
write.csv(Reduce("+", corr.data.X.rf) / length(corr.data.X.rf), file = "Simulation_13_ref_dat_corr.csv")
