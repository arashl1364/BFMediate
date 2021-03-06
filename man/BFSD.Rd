% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BFSD.R
\name{BFSD}
\alias{BFSD}
\title{log Bayes factor for full mediation}
\usage{
BFSD(Post,Prior,burnin)
}
\arguments{
\item{Post}{output from PartialMed or any of measurement models}

\item{Prior}{prior variance of the direct effect}

\item{burnin}{number of MCMC draws before the posterior is converged, default = R/5}
}
\value{
log(BF_01), which is the evidence in favor of/against full mediation model (see Laghaie and Otter (2020) for guidelines on how to interpret BF_01)
}
\description{
Computes log Bayes factors for the partial mediation model using Savage-Dickey ratio
}
\examples{
simPartialMed = function(beta_M,beta_Y, sigma_M, sigma_Y,N,X) {
  eps_M = rnorm(N)*sigma_M # generate errors for M (independent)
  eps_Y = rnorm(N)*sigma_Y # generate errors for Y (independent)
  M = beta_M[1] + beta_M[2] * X + eps_M # generate latent mediator M
  Y = beta_Y[1] + beta_Y[2] * M + beta_Y[3] * X + eps_Y # generate dependent variable
  list(X = X, M = M, Y = Y)
}

# Set up data generating parameters
N = 1000 # number of observations
sigma_M = .2^.5 # error std M
sigma_Y = .2^.5 # error std Y
beta_M = c(1, .3) # beta_0M and beta_1
beta_Y = c(1, .5, 0) # beta_0Y, beta_2, beta_3
X = rnorm(N,mean = 1,sd = 1)# generate random X
# Generate data based on parameters
Data = simPartialMed(beta_M,beta_Y,sigma_M,sigma_Y,N,X)

#Estimation
A_M = c(100,100); #Prior variance for beta_0M, beta_1
A_Y = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3
R = 2000
out = PartialMed(Data=Data, Prior = list(A_M=A_M, A_Y=A_Y), R = R)
#Computing Bayes factor
BFPartialMed = exp(BFSD(Post = out , Prior = A_Y[3], burnin = R/5))
}
\seealso{
For simulating data from simple mediation model see \link[BFMediate]{PartialMed}
}
