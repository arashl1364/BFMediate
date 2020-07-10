#' Bayes factor for full mediation
#'
#' @description
#' Computes Bayes factors for the partial mediation model using Savage-Dickey approximation
#'
#' @usage BFSD(post,prior,burnin)
#'
#' @param post  output from PartialMed or any of measurement models
#' @param prior  prior variance of the direct effect
#' @param burnin number of MCMC draws before the posterior is converged, default = R/5
#'
#' @return
#' log(BF_01), which is the evidence in favor of the full mediation model (see Laghaie and Otter (2020) for guidelines on how to interpret BF_01)
#' @export
#' @examples
#' # Estimation
#' A_M = c(100,100);    # Prior variance for beta_0M, beta_1
#' A_Y = c(100,100,1)   # Prior variance for beta_0Y, beta_2, beta_3
#' R = 2000
#' out = PartialMed(Data=Data, pars = list(A_M=A_M, A_Y=A_Y), R = R)
#' #Computing Bayes factor
#' BFPartialMed = exp(BFSD(post = out , prior = A_Y[3], burnin = R/5))
#' @seealso
#' For simulating data from simple mediation model see \link[BFMediate]{PartialMed}
#Description:
# BFSD computes Bayes factors for the partial mediation model using Savage-Dickey approximation.
# The restricted model is full mediation (direct effect = 0) and the unrestricted model (direct effect != 0) )
# is partial mediation.
#Arguments:
# post  output from PartialMed or any of measurement models
# prior prior variance of the direct effect
# burnin number of MCMC draws before the posterior is converged (def: R/5 )
#Details:
#
#Value:
# log(BF_01), which is the evidence in favor of the full mediation model (see Laghaie and Otter (2020) for guidelines on how to interpret BF_01)
BFSD = function(post,prior=1,burnin){ #function(Data,post,prior,model,burnin,M){

  A = prior
  lik=0;

  ######################################
  ######################################
  R =  length(post$mu_draw)        #length(post$ssq_M)
  # k=dim(as.matrix(post$beta_1))[2] #for stan computations k is =1 in unidimensional X case, and >1 for multidimensional x, for runiregGibbs (since we estimate intercept, k will be =2 and higher)
  # Data$X = cbind(rep(1,nobs),Data$X)
  outer= rep(0,(R-burnin)) #the outer sum

  #computing the marginal posterior of beta_3
  for(r in (burnin+1):R){
    outer[r-burnin] = dnorm(x = 0, mean = post$mu_draw[r], sd = sqrt(post$var_draw[r]), log = T)
    # outer[r-burnin] = dnorm(x = 0, mean = post$mubeta_2_draw[r,k+1], sd = sqrt(post$varbeta_2_draw[k+1,k+1,r]), log = T)
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
