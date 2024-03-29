# Mediation Analysis using Latent Variable Models, Bayesian Estimation, and Bayes Factors 

The focus is on measuring evidence for (full) mediation using Bayes factors. The package covers models with measurement error and discretization in the mediator (M) and/or the dependent variable. 

This is a preliminary version of the package.

## Installation 

```
## Installation from github might take a while on some machines
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
devtools::install_github("arashl1364/BFMediate")
```

#  Measuring data based evidence of mediation 

The following examples simulate data, call package functions, and summarize output to illustrate the models and analyses discussed in detail in Laghaie & Otter (2023).


## Simple mediation model 

This example simulates data from a simple mediation model in which the data generating direct effect is zero, i.e., a model with full mediation. 
```
library(BFMediate)

simPartialMed = function(beta_1,beta_2, sigma_M, sigma_Y,N,X) {
  eps_M = rnorm(N)*sigma_M
  eps_Y = rnorm(N)*sigma_Y
  M = beta_1[1] + beta_1[2] * X + eps_M # generate mediator M
  Y = beta_2[1] + beta_2[2] * M + beta_2[3] * X + eps_Y # generate dependent variable
  list(X = X, M = M, Y = Y)
}

N = 1000 # number of observations
sigma_M = 1^.5 # error std M
sigma_Y = 1^.5 # error std Y
beta_1 = c(1, .5) # beta_0M and beta_1
beta_2 = c(1, 1.5, 0) # beta_0Y, beta_2, beta_3
X = rnorm(N,mean = 1,sd = 1) # generate random X
# generate data based on parameters
Data = simPartialMed(beta_1,beta_2,sigma_M,sigma_Y,N,X)

```
We next fit a simple mediation model to the data just generated and compute the Bayes factor for the model, using the function _Mediate_. The function takes Data (here, a list containing the manipulation (X), the mediator (M), and the dependent variable (Y)), Prior (a list of two vectors of prior variances in the first and the second mediation equation, A_M and A_Y, organized in list form), the number of Mcmc iterations R, and burnin (the number of initial Mcmc draws that are excluded when computing the Bayes factor). 

```
# Estimation and Bayes factor computation

# Choosing the reference prior for the direct effect (beta_3)
A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y_ref = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3(reference prior)
Prior = list(A_M=A_M, A_Y=A_Y_ref)
R = 10000
burnin = 2000
out = Mediate(Data = Data, Model = "Simple",Prior = Prior,R = R, burnin = burnin)

```

The _Mediate_ function facilitates the comparison of mediation analysis results across different measurement approaches, i.e., Baron & Kenny (1986), Preacher & Hayes (2004), and the proposed Bayesian approach to mediation analysis.  

```

# Regression estimates of B&K's M and Y equations
out$BK$eq1
out$BK$eq2
out$BK$FullMed    # full mediation Null hypothesis test result

# Bootstrapped estimates a la Preacher & Hayes (2004)
out$PH$Indirect_mean  # indirect effect's mean
out$PH$Indirect_CI    # indirect effect's bootstrapped confidence interval
out$PH$Direct_CI      # direct effect's bootstrapped confidence interval

# Proposed Bayesian approach results
out$Simple$Indirect_CI    # indirect effect 95% posterior credible interval
out$Simple$Direct_CI      # direct effect 95% posterior credible interval
out$Simple$BF             # Bayes factor
out$Simple$evidence       # evidence in favor of full mediation (Kass & Raftery 1995)

```
The NHST test result from the Baron & Kenny (1986) procedure fails to reject full mediation, and the bootstrapped direct effect includes zero. However, as dicussed in Laghaie & Otter (2023), these could result from low power and high sampling variability, and are not necessarily evidence for mediation. Computing Bayes factor on the other hand provide a measure of data based evidence supporting/against mediation. The output of _Mediate_ also contains all the parameter estimates (posterior draws). For a complete list of the output values as well as the prior specification for all the model parameters please see the help page of the function.


## Model with multiple (discretized) indicators for M and Y

If, as often the case, multiple indicators for M and Y are available in the form of scale ratings, they can be used to control for both measurement error and discretization in a latent variable model(for a visual illustration of the model please see the web appendix of Laghaie & Otter (2023)). Based on the latent variable model, Bayes factors then test for conditional mean independence between X and (latent) Y given (latent) M.  Note that conditional mean independence at the latent variable level may (strongly) hold even if conditional mean independence is (strongly) rejected at the level of observed variables. We illustrate this result with simulated data from a full mediation model with multiple categorized indicators:      

```
set.seed(60)
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
                    -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Mcut, byrow = T)  #cutoffs for discretization of M indicators

cutoff_Y =  matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
                     -100, 0, 1, 2, 3, 4, 5, 100),ncol= Ycut, byrow = T)  #cutoffs for discretization of Y indicators
DataMYCat = SimMeasurementMYCat(X, beta_1, cutoff_M, beta_2, cutoff_Y, M_ind, Y_ind, lambda, tau, ssq_m_star, ssq_y_star)

```

## Estimating a simple model with composite measures 

We first average out the indicators to make composite measures and estimate a simple model using the those measures. 

```
#creating composite measures of M and Y using their indicators
Data_comp = NULL
Data_comp$M = rowMeans(DataMYCat$m_tilde)
Data_comp$Y = rowMeans(DataMYCat$y_tilde)
Data_comp$X = DataMYCat$X

#Estimation
A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3(reference prior)
Prior = list(A_Y = A_Y, A_M = A_M)
out_comp = Mediate(Data = Data_comp, Model = "Simple",Prior = Prior,R = R, burnin = burnin)

```

We compare the simple model analysis results across different procedures. 

```
#Regression estimates of B&K's M and Y equations
out_comp$BK$eq1
out_comp$BK$eq2
out_comp$BK$FullMed #full mediation Null hypothesis test result

#Bootstrapped estimates a la Preacher & Hayes (2004)
out_comp$PH$Indirect_mean  #indirect effect's mean
out_comp$PH$Indirect_CI    #indirect effect's bootstrapped confidence interval
out_comp$PH$Direct_CI     #direct effect's bootstrapped confidence interval

#Proposed Bayesian approach results
colMeans(out_comp$Simple$beta_M)   #Bayesian posterior means of M equation's parameters
colMeans(out_comp$Simple$beta_Y)   #Bayesian posterior means of Y equation's parameters
out_comp$Simple$BF     #Bayes factor
out_comp$Simple$evidence    #evidence in favor of full mediation (Kass & Raftery 1995)

```

Baron & Kenny (1986) and bootstrapping methods both reject full mediation and even though the indirect effect is significant, we cannot rule out alternative models. Also the Bayes factor of the simple model shows strong evidence against full mediation. 

##  Estimating a latent variable model

Instead of estimating a simple model using the composite measures, we can estimate the LVM to account for measurement error and descritization using _Mediate_ by setting Model = "MYCat". Note that Data argument here should be a list containing 
- X a vector of the manipulated variable; 
- m_tilde and y_tilde, matrices of M and Y indicator variables respectively, stored column-wise.

```
out_MYCat = Mediate(Data = DataMYCat, Model = "MYCat",Prior = Prior,R = 10000, burnin = 2000)

```

The comparison of results from the LVM to those from the simple model illustrates how accounting for measurement error improves inference for the underlying model structure and the direct and indirect effect:

```
out_MYCat$Indirect_CI    
out_MYCat$Direct_CI      
out_MYCat$BF
out_MYCat$evidence

```

# BF sensitivity analysis and MCMC convergence

## BF sensitivity to the choice of prior and convergence in the simple mediation model

For valid inference, e.g., about posterior means, 95% CIs, and meaningful Bayes factors computed from the posterior draws, it is important ensure to the MCMC has converged.  One way to inspect convergence is to plot MCMC traces, i.e., the MCMC-iteration count on the x-axis agains draws of functions of draws on the y-axis.
To illustrate the convergence we simulates data from a model with full mediation. 

```

N = 1000 # number of observations
sigma_M = 1^.5 # error std M
sigma_Y = 1^.5 # error std Y
beta_1 = c(1, .5) # beta_0M and beta_1
beta_2 = c(1, 1.5, 0) # beta_0Y, beta_2, beta_3
X = rnorm(N,mean = 1,sd = 1) # generate random X
# generate data based on parameters
Data = simPartialMed(beta_1,beta_2,sigma_M,sigma_Y,N,X)

```
We next estimate a simple mediation model conditional on the data just generated, using the function _PartialMed_.

```
#Estimation

#Choosing the reference prior for the direct effect (beta_3)
A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y_ref = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3(reference prior)
R = 2000
out_1 = PartialMed(Data=Data, Prior = list(A_M=A_M, A_Y=A_Y_ref), R = R)

```
Even though the simple mediation model employs a Gibbs sampler, it essentially converges immediately. More specifically, based on R=2000 no burn-in period from starting values can be discerned:

```
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}

#plotting MCMC traces to assess convergence
quartz()
par(mfrow=c(1,3)) 
matplot(out_1$beta_M[,2], type = 'l')
matplot(out_1$beta_Y[,2], type = 'l')
matplot(out_1$beta_Y[,3], type = 'l')
```
As an example of a posterior summary of interest posterior mean estimates of regression coefficients are computed next:

```
#estimation results
colMeans(out_1$beta_M)    #posterior means of M equation's coeffcients
colMeans(out_1$beta_Y)    #posterior means of Y equation's coeffcients

```

The distribution of the indirect effect is obtained by multiplying posterior draws of beta_1 and beta_2. An assessment of marginal posterior uncertainty is available by computing estimates of posterior quantiles:    

```
hist(out_1$beta_M[,2] * out_1$beta_Y[,2])   #plotting the distribution of the indirect effect
quantile(out_1$beta_M[,2]*out_1$beta_Y[,2], probs = c(.025,.975))  #95% posterior credible interval(CI) of the indirect effect
```

Next we re-estimate a model with a much more diffuse prior for the direct effect (the thirs element of A_Y is set to 100 below) to illustrate the sensitivity of the Bayes factor testing conditional mean independence to 
this subjective prior choice:

```
#Choosing a diffuse prior for the direct effect (beta_3)
A_Y_dif = c(100,100,100) #Prior variance for beta_0Y, beta_2, beta_3
R = 2000
out_100 = PartialMed(Data=Data, Prior = list(A_M=A_M, A_Y=A_Y_dif), R = R)

```
The function _BFSD_ computes log Bayes factors for each model.  The resulting Bayes factors measure the degree of empirical support for the Null hypothesis that the direct effect is equal to zero, i.e., that X and Y are mean independent conditional on M. For the prior specification of the rest of the model parameters please see the help page of the _PartialMed_ function.

```
#comparing  Bayes factors of the two models
BF_1 = exp(BFSD(Post = out_1 , Prior = A_Y_ref[3], burnin = 0))
BF_100 = exp(BFSD(Post = out_100 , Prior = A_Y_dif[3], burnin = 0))

```

The Bayes factor of the model with the diffuse prior (prior variance equal to 100) is much larger than the one that uses the suggested reference prior (prior variance equalt to 1).  This is because a diffuse prior places a relatively higher probability on (absolutely) very large direct effects.


# Sensitivity of variance explained to the choice of prior variance and convergence in the model with multiple (discretized) indicators for M and Y 

The model with discretized ordinal indicators is identified by fixing the residual variances of
latent M and Y as well as slope parameters connecting latent M and Y to latent continuous
indicators to 1. Under this identification strategy, the natural target for assessing the
subjective prior for the direct effects from X to latent Y is the prior (conditional) variance
explained in latent Y. 
To illustrate the prior assessment and convergence we simulate data from a full mediation model with multiple categorized indicators:      

```

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
                    -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Mcut, byrow = T)  #cutoffs for discretization of M indicators

cutoff_Y =  matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
                     -100, 0, 1, 2, 3, 4, 5, 100),ncol= Ycut, byrow = T)  #cutoffs for discretization of Y indicators
DataMYCat = SimMeasurementMYCat(X, beta_1, cutoff_M, beta_2, cutoff_Y, M_ind, Y_ind, lambda, tau, ssq_m_star, ssq_y_star)

```

We assess the subjective prior for the direct effect prior by plotting the corresponding distribution of variance explained in the dependent variable: 

```
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
```
The histogram shows that the reference prior results in a distribtution of variance explained concentrated away from 1, i.e., embodies the idea that direct effects are unlikely to dominate the model a priori.

The estimation function _MeasurementMYCat_ takes the following arguments: 
Data: a list containing 
- X a 2-column matrix with a vector of ones as its first and the manipulated variable as its second column; 
- m_tilde and y_tilde, matrices of M and Y indicator variables respectively, stored column-wise; - k_M and k_Y, scalars, the number of rating scale points for M and Y indicators respectively; 
- M_ind and Y_ind, scalars, the number of M and Y indicators.
Prior: a list containing vectors of prior variances for regression coefficients in the first and second (latent) mediation equations.  
R: number of Mcmc iterations.
```
Mcut = max(DataMYCat$m_tilde) + 1
Ycut = max(DataMYCat$y_tilde) + 1
Data = list(X=cbind(rep(1,length(DataMYCat$X)),DataMYCat$X), m_tilde=as.matrix(DataMYCat$m_tilde),
            y_tilde=as.matrix(DataMYCat$y_tilde), k_M = Mcut-1, k_Y=Ycut-1,
            M_ind=dim(as.matrix(DataMYCat$m_tilde))[2], Y_ind=dim(as.matrix(DataMYCat$y_tilde))[2])

A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3(reference prior)

```
## Estimation and convergence assessment 
To illustrate convergence of the MCMC, and the assessment of convergence, we estimate the model with different numbers of iterations (R), and plot selected posterior MCMC traces.

```
R = 2000
out_2000 = MeasurementMYCat(Data=Data, Prior = list(A_M=A_M, A_Y=A_Y), R = R)
#plotting MCMC traces to assess convergence
quartz()
par(mfrow=c(1,3)) 
matplot(out_2000$beta_M[,2], type = 'l')
matplot(out_2000$beta_Y[,2], type = 'l')
matplot(out_2000$beta_Y[,3], type = 'l')

R = 10000
out_10000 = MeasurementMYCat(Data=Data, Prior = list(A_M=A_M, A_Y=A_Y), R = R)
#plotting MCMC traces to assess convergence
quartz()
par(mfrow=c(1,3)) 
matplot(out_10000$beta_M[,2], type = 'l')
matplot(out_10000$beta_Y[,2], type = 'l')
matplot(out_10000$beta_Y[,3], type = 'l')

```
Based in R = 2000, we see that posteriors appear to converge quite fast (after 500 iterations or so). To ensure precise posterior summaries, we estimate the model with R = 10000 and throw away the first 2000 draws as burnin. Run time is not much of a constraint here.


Based on on 8000 posterior draws, one can now compute any function of the parameters and summarize its distribution, such as e.g., the direct and indirect effects.

```
#results
burnin = 2000
colMeans(out_10000$beta_M[(burnin+1):R,])    #posterior means of M equation's coeffcients
colMeans(out_10000$beta_Y[(burnin+1):R,])    #posterior means of Y equation's coeffcients
apply(out_10000$lambdadraw[,,(burnin+1):R],MARGIN = c(1,2),FUN = mean) #posterior means of the M measurement equation parameters (column1: intercepts, column2: coefficients(fixed to 1))
apply(out_10000$taudraw[,,(burnin+1):R],MARGIN = c(1,2),FUN = mean)    #posterior means of the Y measurement equation parameters (column1: intercepts, column2: coefficients(fixed to 1))

apply(out_10000$cutoff_M[,,(burnin+1):R],c(1,2),FUN = mean)   #posterior means of M indicators' cutoffs
apply(out_10000$cutoff_Y[,,(burnin+1):R],c(1,2),FUN = mean)   #posterior means of Y indicators' cutoffs
```


