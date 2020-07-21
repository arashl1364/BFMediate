library(devtools)
install_github("arashl1364/BFMediate")
library(BFMediate)
##################################################
######   Simple partial mediation model    #######
#### and BF sensitivity to the choice of prior ###
##################################################
simPartialMed = function(beta_1,beta_2, sigma_M, sigma_Y,N,X) {
  eps_M = rnorm(N)*sigma_M
  eps_Y = rnorm(N)*sigma_Y
  M = beta_1[1] + beta_1[2] * X + eps_M # generate mediator M
  Y = beta_2[1] + beta_2[2] * M + beta_2[3] * X + eps_Y # generate dependent variable
  list(X = X, M = M, Y = Y)
}

#Data generation
N = 1000 # number of observations
sigma_M = 1^.5 # error std M
sigma_Y = 1^.5 # error std Y
beta_1 = c(1, .5) # beta_0M and beta_1
beta_2 = c(1, 1.5, 0) # beta_0Y, beta_2, beta_3
X = rnorm(N,mean = 1,sd = 1) # generate random X
# generate data based on parameters
Data = simPartialMed(beta_1,beta_2,sigma_M,sigma_Y,N,X)

#Estimation

#Choosing the reference prior for the direct effect (beta_3)
A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y_inf = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3(reference prior)
R = 2000
out_1 = PartialMed(Data=Data, Pars = list(A_M=A_M, A_Y=A_Y_inf), R = R)
#estimation results
colMeans(out_1$beta_1)    #posterior means of the first mediation equation's coeffcients
colMeans(out_1$beta_2)    #posterior means of the second mediation equation's coeffcients


#Choosing a diffuse prior for the direct effect (beta_3)
A_Y_dif = c(100,100,100) #Prior variance for beta_0Y, beta_2, beta_3
R = 2000
out_100 = PartialMed(Data=Data, Pars = list(A_M=A_M, A_Y=A_Y_dif), R = R)
#estimation results
colMeans(out_1$beta_1)    #posterior means of the first mediation equation's coeffcients
colMeans(out_1$beta_2)    #posterior means of the second mediation equation's coeffcients


#comparing  Bayes factors of the two models
BF_1 = exp(BFSD(Post = out_1 , Prior = A_Y_inf[3], burnin = 0))
BF_100 = exp(BFSD(Post = out_100 , Prior = A_Y_dif[3], burnin = 0))
#BF_100 is much larger than BF_1 as the former model assumes large values of direct effect
#a priori probable

######################################################################
#####  Model with multiple (continuous) indicators for M and Y    ####
##### and its BF sensitivity to the choice of reference indicator ####
######################################################################
library(rstan)
SimMeasurementCont = function( beta_1, beta_2 , lambda, tau, m_ind, y_ind, sigma_M, sigma_m_star, sigma_y, sigma_y_star, N, X) {

  m_star = matrix(double(N*m_ind),ncol = m_ind); y_star = matrix(double(N*y_ind),ncol = y_ind)
  eps_m_star = matrix(double(N*m_ind),ncol = m_ind); eps_y_star = matrix(double(N*y_ind),ncol = y_ind)
  eps_M = rnorm(N)*sigma_M # generate errors for M (independent)
  eps_Y = rnorm(N)*sigma_y # generate errors for y (independent)
  M = beta_1[1] + X*beta_1[2]  + eps_M # generate latent mediator M
  y = beta_2[1] +  M*beta_2[2] + X*beta_2[3] + eps_Y # generate dependent variable

  eps_m_star[,1]=rnorm(N)*sigma_m_star[1] # generate errors for m_star (independent)
  m_star[,1] =  M + eps_m_star[,1] # generate observed mediator indicators m_star
  if(m_ind>1){
    for(i in 2:(m_ind))   {
      eps_m_star[,i]=rnorm(N)*sigma_m_star[i] # generate errors for m_star (independent)
      m_star[,i] =  lambda[(i-1),1] + M*lambda[(i-1),2]  + eps_m_star[,i]
    }
  }

  eps_y_star[,1]=rnorm(N)*sigma_y_star[1] # generate errors for y_star (independent)
  y_star[,1] =  y + eps_y_star[,1] # generate observed dependent variable indicators y_star
  if(y_ind>1){
    for(i in 2:(y_ind)){
      eps_y_star[,i]=rnorm(N)*sigma_y_star[i] # generate errors for y_star (independent)
      y_star[,i] =  tau[(i-1),1] + y*tau[(i-1),2] + eps_y_star[,i]
    }
  }
  list(X = X, M = M, m_star = m_star, y = y, y_star=y_star)
}

#Data generation
m_ind = 2; y_ind = 3;
sigma_M = 1^.5 # error std M
sigma_y = 1^.5 # error std y
sigma_m_star = c(1,.5)^.5 #c(1,2)^.5
sigma_y_star = c(1,.5,.5)^.5  #c(2,1)^.5
beta_1 = c(1,1)
beta_2 = c(1,3,0)
lambda = matrix(c(.5,1),ncol=2,byrow = T)  #tau and lambda should have y_ind-1 and m_ind-1 rows and 2 columns
tau =  matrix(c(.5,.1,
                .5,.3),ncol=2,byrow = T)
k=length(beta_1)-1
nobs = 500   # number of observations
X = runif(nobs) # generate random X from a uniform distribution
Data = SimMeasurementCont( beta_1, beta_2 , lambda, tau, m_ind, y_ind, sigma_M, sigma_m_star, sigma_y, sigma_y_star, nobs, X)

R = 5000; burnin = 3000
A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3(reference prior)

#Estimation with the largest y indicator as the identifying indicator
out_large_ind = MeasurementCont(Data = Data, Prior = list(A_M = A_M, A_Y = A_Y),R=5000, burnin = 3000)

#results
colMeans(out_large_ind$beta_1)    #posterior means of the first mediation equation's coeffcients
colMeans(out_large_ind$beta_2)    #posterior means of the second mediation equation's coeffcients
colMeans(out_large_ind$tau)     #posterior means of the M measurement equation parameters (column1: intercepts, column2: coefficients(row1 is fixed to {0,1}))
colMeans(out_large_ind$lambda)  #posterior means of the M measurement equation parameters (column1: intercepts, column2: coefficients(row1 is fixed to {0,1}))

#Estimation with the smallest y indicator as the identifying indicator
temp = Data$y_star[,1]; Data$y_star[,1] = Data$y_star[,2]; Data$y_star[,2] = temp   #changing the refrence indicator of y
out_small_ind = MeasurementCont(Data = Data, Prior = list(A_M = A_M, A_Y = A_Y),R=5000, burnin = 3000)

#Comparing Bayes factors of the two models
BF_large_ind = exp(BFSD(Post = out_large_ind , Prior = A_Y[3],burnin = 0))
BF_small_ind = exp(BFSD(Post = out_small_ind , Prior = A_Y[3],burnin = 0))

####################################################################################
########## Model with multiple (categorical) indicators for M and Y,    ############
##### and the sensitivity of variance explained to the choice of prior variance  ###
####################################################################################
SimMeasurementMYCat = function(X, beta_1, cutoff_M, beta_2, cutoff_Y, M_ind, Y_ind, lambda, tau, ssq_m_star, ssq_y_star){

  nobs = dim(X)[1]
  m_tilde = m_star = matrix(double(nobs*M_ind), ncol = M_ind)
  y_tilde = y_star = matrix(double(nobs*Y_ind), ncol = Y_ind)

  M = cbind(rep(1,nobs),X)%*%beta_1 + rnorm(nobs)
  for(i in 1: M_ind){
    m_star[,i] = lambda[i] + M + sqrt(ssq_m_star[i])*rnorm(nobs);
    m_tilde[,i] = cut(m_star[,i], br = cutoff_M[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
  }

  Y = cbind(rep(1,nobs),cbind(M,X))%*%beta_2 + rnorm(nobs)
  for(i in 1: Y_ind){
    y_star[,i] = tau[i] + Y + sqrt(ssq_y_star[i])*rnorm(nobs);
    y_tilde[,i] = cut(y_star[,i], br = cutoff_Y[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
  }

  return(list(Y = Y, M = M, y_tilde = y_tilde, m_tilde = m_tilde, X = X,
              k_M=dim(cutoff_M)[2]-1, beta_1 = beta_1, lambda = lambda, ssq_m_star = ssq_m_star, m_star = m_star, cutoff_M = cutoff_M,
              k_Y=dim(cutoff_Y)[2]-1, beta_2 = beta_2, tau = tau, ssq_y_star = ssq_y_star, y_star = y_star, cutoff_Y = cutoff_Y,
              M_ind=M_ind, Y_ind=Y_ind))
}


#Data generation
M_ind = 2   #number of M indicators
Y_ind = 2   #number of Y indicators
Mcut = Ycut = 8  #number of cutoffs (including the -inf and inf which are set -100 and 100 here)
nobs=1000    #number of observations
X=as.matrix(runif(nobs,min=0, max=1))  #generating X
beta_1 = c(1,.8)
beta_2 = c(1,1.5, 0)
ssq_m_star = c(.5,.3)   #measurement error variances of M indicators
lambda = c(0,-.5)  #the intercepts of the latent M indicators w. measurement error
ssq_y_star = c(.2,.2)  #measurement error variances of Y indicators
tau = c(0,-.5)   #the intercepts of the latent Y indicators w. measurement error
cutoff_M = matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
                    -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Mcut, byrow = T)

cutoff_Y =  matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
                     -100, 0, 1, 2, 3, 4, 5, 100),ncol= Ycut, byrow = T)
DataMYCat = SimMeasurementMYCat(X, beta_1, cutoff_M, beta_2, cutoff_Y, M_ind, Y_ind, lambda, tau, ssq_m_star, ssq_y_star)

##Sensitivity analysis of variance explained to the choice of prior
sd = c(sqrt(.1),1,sqrt(10),10)
var_X = var(X)
ssq_y = 1
quartz()
par(mfrow=c(2,2))
var_explained = matrix(double(1000*4), ncol = 4)
for(i in 1:length(sd)){

  beta_3_prior = rnorm(1000,sd = sd[i])
  direct_var1 = (beta_3_prior^2) * c(var_X)
  direct_var = direct_var1
  indirect_var = beta_2[2]^2 * (beta_1[2]^2 * c(var(X)) + 1)
  var_explained[,i] = direct_var/(direct_var+indirect_var+ssq_y)

  hist(var_explained[,i],main= paste("Var = ",sd[i]^2))
}

#Estimation
Mcut = max(DataMYCat$m_tilde) + 1
Ycut = max(DataMYCat$y_tilde) + 1
Data = list(X=cbind(rep(1,length(DataMYCat$X)),DataMYCat$X), m_tilde=as.matrix(DataMYCat$m_tilde),
            y_tilde=as.matrix(DataMYCat$y_tilde), k_M = Mcut-1, k_Y=Ycut-1,
            M_ind=dim(as.matrix(DataMYCat$m_tilde))[2], Y_ind=dim(as.matrix(DataMYCat$y_tilde))[2])

A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3(reference prior)

out = MeasurementMYCat(Data=Data, Prior = list(A_M=A_M, A_Y=A_Y), R = 10000)

#results
colMeans(out$beta_1)    #posterior means of the first mediation equation's coeffcients
colMeans(out$beta_2)    #posterior means of the second mediation equation's coeffcients
apply(out$lambdadraw,MARGIN = c(1,2),FUN = mean) #posterior means of the M measurement equation parameters (column1: intercepts, column2: coefficients(fixed to 1))
apply(out$taudraw,MARGIN = c(1,2),FUN = mean)    #posterior means of the Y measurement equation parameters (column1: intercepts, column2: coefficients(fixed to 1))

apply(out$cutoff_M,c(1,2),FUN = mean)   #posterior means of M indicators' cutoffs
apply(out$cutoff_Y,c(1,2),FUN = mean)   #posterior means of Y indicators' cutoffs


