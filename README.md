# Mediation Analysis using Latent Variable Models and Bayesian Estimation, and Bayes factors 

The focus is on measuring evidence for (full) mediation using Bayes factors. The package covers models with measurement error and discretization in the mediator (M) and/or the dependent variable. 

Note: This is a preliminary version of the package. Please do not quote or distribute.

## Installation 

```
## Installation from github might take a while on some machines
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
devtools::install_github("arashl1364/BFMediate")
```

## Note: Compiling the Stan part of the package might require the following configuration on some Windows machines. If the installation fails please run the code below and try again: 

```
dotR <- file.path(Sys.getenv("HOME"), ".R")
if (!file.exists(dotR)) dir.create(dotR)
M <- file.path(dotR, "Makevars.win")
if (!file.exists(M)) file.create(M)
cat("\nCXX14FLAGS=-O3 -march=corei7 -mtune=corei7",
    "CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
    "CXX11FLAGS=-O3 -march=corei7 -mtune=corei7",
    file = M, sep = "\n", append = TRUE)
```

## Getting Started

The following examples simulate data, call package functions, and summarize output to illustrate the models and analyses discussed in detail in Laghaie & Otter (2020).

# Simple partial mediation model and BF sensitivity to the choice of prior

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
We next estimate a simple mediation model conditional on the data just generated, using the function _PartialMed_.  The function takes Data (a list containing the manipulation (X), the mediator (M), and the dependent variable (Y)), Prior (two vectors of prior variances in the first and the second mediation equation, A_M and A_Y, organized in list form), and the number of Mcmc iterations R. 

```
#Estimation

#Choosing the reference prior for the direct effect (beta_3)
A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y_ref = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3(reference prior)
R = 2000
out_1 = PartialMed(Data=Data, Prior = list(A_M=A_M, A_Y=A_Y_ref), R = R)

```
For valid inference, e.g., about posterior means, 95% HDIs, and meaningful Bayes factors computed from the posterior draws, it is important ensure the MCMC has converged.  One way to inspect convergence is to plot MCMC traces, i.e., the MCMC-iteration count on the x-axis agains draws of functions of draws on the y-axis.
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

Distribution of the indirect effect can be computed simply by multiplying the posteriors of beta_1 and beta_2. A simple assessment of marginal posterior uncertainty is available by computing estimates of posterior quantiles:   

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

# Model with multiple (discretized) indicators for M and Y, and the sensitivity of variance explained to the choice of prior variance  

If, as often the case, multiple indicators for M and Y are available in the form of scale ratings, the function _MeasurementMYCat_ uses these indicators to control for both measurement error and discretization in a latent variable model. Based on the latent variable model, Bayes factors then test for conditional mean independence between X and (latent) Y given (latent) M.  Note that conditonal mean independence at the latent variable level may (strongly) hold even if conditional mean independence is (strongly) rejected at the level of observed variables.  We illustrate this result with simulated data from a full mediation model with multiple categorized indicators:      

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

We assess the subjective prior for the direct effect prior by plotting the corresponding distribution of  variance explained in the dependent variable: 

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

## Comparing results from different approaches 

The _Mediate_ function in the package facilitates the comparison of mediation analysis results across different measurement approaches, i.e., Baron & Kenny (1986), Preacher & Hayes (2004), and the proposed Bayesian approach to mediation analysis. Next is an illustration using the data from the previous example.

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

If "Simple" is given as the model, _Mediate_ function performs a mediation analysis a la Baron & Kenny (1986), as well as, Preacher & Hayes (2004) bootstrapping analysis and the proposed Bayesian method. The results are stored in lists BK, PH, and Simple respectively.

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

Instead of estimating a simple model using the composite measures, we can estimate the LVM to account for measurement error and descritization by setting Model = "MYCat". Note that here the analysis done by _Mediate_ essentially runs _MeasurementMYCat_ and _BFSD_, as illustrated previously. 

```
out_MYCat = Mediate(Data = DataMYCat, Model = "MYCat",Prior = Prior,R = 10000, burnin = 2000)

```

The comparison of results from the LVM to those from the simple model illustrates how accounting for measurement error improves inference for the underlying model structure and the direct and indirect effect:

```
out_MYCat$BF
out_MYCat$evidence
```
Draws for all the parameters are stored in the output. For a complete list of the output values as well as the prior specification for all the model parameters please see the help page of the function.   

## Stability of Bayes factor 
Last, we illustrate how to assess Monte-Carlo variation in inferred Bayes factors.  Based on perfect knowledge of the posterior the Bayes factor is a fixed number given the model and the prior.  Because inferene is based on a finite sample from the posterior, it is important to check if the posterior sample is large enough to yield the same Monte-Carlo estimate of the Bayes factor (up to a required degree of precision) across independent replications.  

```
rep = c(10,50,100,1000,2000,5000,8000)
BF = rep(0,length(rep))
samp = NULL
for(i in 1:length(rep)){
  draws = sample((burnin+1):R, size = rep[i])
  samp$mu_draw = out_MYCat$mu_draw[draws]
  samp$var_draw = out_MYCat$var_draw[draws]
  BF[i] = exp(BFSD(samp,Prior = 1, burnin = 0))
}
BF
```

The above analysis shows that Monte-Carlo based estimates of Bayes stabilize already with posterior samples of size 2000.  However, this is a somewhat optimistic estimate because it is based on randomly drawn 2000 posterior draws from a larger posterior sample.  An alternative assessment is based on batch estimates and their stability across different batch sizes.


