#' Bayes factor for full mediation
#'
#' @description
#' Computes Bayes factors for the partial mediation model using Savage-Dickey approximation
#'
#' @usage BFSD(Post,Prior,burnin)
#'
#' @param Post  output from PartialMed or any of measurement models
#' @param Prior  prior variance of the direct effect
#' @param burnin number of MCMC draws before the posterior is converged, default = R/5
#'
#' @return
#' log(BF_01), which is the evidence in favor of the full mediation model (see Laghaie and Otter (2020) for guidelines on how to interpret BF_01)
#' @export
#' @examples
#' simPartialMed = function(beta_1,beta_2, sigma_M, sigma_Y,N,X) {
#'   eps_M = rnorm(N)*sigma_M # generate errors for M (independent)
#'   eps_Y = rnorm(N)*sigma_Y # generate errors for Y (independent)
#'   M = beta_1[1] + beta_1[2] * X + eps_M # generate latent mediator M
#'   Y = beta_2[1] + beta_2[2] * M + beta_2[3] * X + eps_Y # generate dependent variable
#'   list(X = X, M = M, Y = Y)
#' }
#'
#' # Set up data generating parameters
#' N = 1000 # number of observations
#' sigma_M = .2^.5 # error std M
#' sigma_Y = .2^.5 # error std Y
#' beta_1 = c(1, .3) # beta_0M and beta_1
#' beta_2 = c(1, .5, 0) # beta_0Y, beta_2, beta_3
#' X = rnorm(N,mean = 1,sd = 1)# generate random X
#' # Generate data based on parameters
#' Data = simPartialMed(beta_1,beta_2,sigma_M,sigma_Y,N,X)
#'
#' #Estimation
#' A_M = c(100,100); #Prior variance for beta_0M, beta_1
#' A_Y = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3
#' R = 2000
#' out = PartialMed(Data=Data, pars = list(A_M=A_M, A_Y=A_Y), R = R)
#' #Computing Bayes factor
#' BFPartialMed = exp(BFSD(Post = out , Prior = A_Y[3], burnin = R/5))
#' @seealso
#' For simulating data from simple mediation model see \link[BFMediate]{PartialMed}
#Description:
# BFSD computes Bayes factors for the partial mediation model using Savage-Dickey approximation.
# The restricted model is full mediation (direct effect = 0) and the unrestricted model (direct effect != 0) )
# is partial mediation.
#Arguments:
# Post  output from PartialMed or any of measurement models
# Prior prior variance of the direct effect
# burnin number of MCMC draws before the posterior is converged
#Details:
#
#Value:
# log(BF_01), which is the evidence in favor of the full mediation model (see Laghaie and Otter (2020) for guidelines on how to interpret BF_01)
BFSD = function(Post,Prior=1,burnin){ #function(Data,Post,Prior,model,burnin,M){

  A = Prior
  lik=0;

  ######################################
  ######################################
  R =  length(Post$mu_draw)        #length(Post$ssq_M)
  # k=dim(as.matrix(Post$beta_1))[2] #for stan computations k is =1 in unidimensional X case, and >1 for multidimensional x, for runiregGibbs (since we estimate intercept, k will be =2 and higher)
  # Data$X = cbind(rep(1,nobs),Data$X)
  outer= rep(0,(R-burnin)) #the outer sum

  #computing the marginal posterior of beta_3
  for(r in (burnin+1):R){
    outer[r-burnin] = dnorm(x = 0, mean = Post$mu_draw[r], sd = sqrt(Post$var_draw[r]), log = T)
    # outer[r-burnin] = dnorm(x = 0, mean = Post$mubeta_2_draw[r,k+1], sd = sqrt(Post$varbeta_2_draw[k+1,k+1,r]), log = T)
  }

  #Numerically Stably computing the log of the outer sum
  outermax = max(outer)
  outer = outer[outer != max(outer)] #removes the maximum element from the vector
  lik = (outermax + log(1 + sum(exp(outer - outermax))))  #the outer sum in the log-scale
  lik = -log((R-burnin)) + lik

  lik = lik - dnorm(x = 0, mean = 0, sd = sqrt(A), log = T)
  lik = -lik #only to adjust to the function output
  ######################################
  ######################################

  return(-lik)
}
