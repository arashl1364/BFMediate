#' Sampler for Partial Mediation Model with Multiple Continuous Indicator for the Mediator and/or DV
#' @description Estimates a partial mediation model with multiple categorical indicators for the mediator and the dependent variable using Hamiltonian Monte Carlo (HMC) with Stan
#'
#' @usage MeasurementCont(Data, Prior, R, burnin)
#'
#' @param Data list(X, m_tilde, y_tilde)
#' @param Prior list(A_M,A_Y)
#' @param R number of MCMC iterations, default = 10000
#'
#' @details
#' ## Model
#'
#' (eq.1) \deqn{M = \beta_0 M + X \beta_1 + U_M}{M = \beta_0*M + X*\beta_1 + U_M}
#' (eq.2) \deqn{Y = \beta_0 Y + M \beta_2 + X \beta_3 + U_Y}{Y = \beta_0*Y + M*\beta_2 + X*\beta_3 + U_Y}
#'
#' ## Indicator equations
#' \deqn{m^*_1 = M + U_{m^*_1}}{m*_1 = M + U_{m*_1}}
#' \deqn{m^*_2 = \lambda_{01} + M + U_{m^*_2}}{m*_2 = \lambda_01 + M + U_{m*_2}}
#' ...
#' \deqn{m^*_k = \lambda_{0k-1} + M + U_{m^*_k}}{m*_k = \lambda_0k-1 + M + U_{m*_k}}
#' \deqn{y^*_1 = M + U_{y^*_1}}{y*_1 = M + U_y*_1}
#' \deqn{y^*_2 = \tau_{01} + M + U_{y^*_2}}{y*_2 = \tau_01 + M + U_{y*_2}}
#' ...
#' \deqn{y^*_l = \tau_{0l-1} + M + U_{y^*_l}}{y*_l = \tau_0l-1 + M + U_{y*_l}}
#'
#'
#' ## Argument Details
#'
#' ## \code{Data = list(X, m_tilde, y_tilde)}
#' \describe{
#'   \item{X(N x 1)}{treatment variable vector}
#'   \item{m_tilde(N x M_ind)}{ediator indicators' matrix }
#'   \item{y_tilde(N x Y_ind)}{dependent variable indicators' matrix}
#' }
#'
#' ## \code{Prior = list(A_M,A_Y)} *\[optional\]*
#' \describe{
#'   \item{A_M}{vector of coefficients' prior variances of eq.1, default = rep(100,2)}
#'   \item{A_Y}{vector of coefficients' prior variances of eq.2, default = c(100,100,1)}
#' }
#'
#' @return
#' \describe{
#'   \item{beta_M(R X 2)}{matrix of eq.1 coefficients' draws}
#'   \item{beta_Y(R X 3)}{matrix of eq.2 coefficients' draws}
#'   \item{lambda (M_ind X 2 X R)}{array of mediator indicator coefficients' draws. Each slice is one draw, where rows represent the indicator equation and columns are the coefficients. All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively.}
#'   \item{ssq_m_star(R X M_ind)}{Matrix of mediator indicator equations' coefficients' error variance draws}
#'   \item{ssq_y_star(R X Y_ind)}{Matrix of dependent variable indicator equations' coefficients' error variance draws}
#'   \item{mu_draw}{vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#'   \item{var_draw}{vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#' }
#'
#' @export
#' @examples
#' library(rstan)
#' SimMeasurementCont = function( beta_1, beta_2 , lambda, tau, m_ind, y_ind, sigma_M, sigma_m_star,
#'                               sigma_y, sigma_y_star, N, X) {
#'
#'   m_star = matrix(double(N*m_ind),ncol = m_ind); y_star = matrix(double(N*y_ind),ncol = y_ind)
#'   eps_m_star = matrix(double(N*m_ind),ncol = m_ind); eps_y_star = matrix(double(N*y_ind),ncol = y_ind)
#'   eps_M = rnorm(N)*sigma_M # generate errors for M (independent)
#'   eps_Y = rnorm(N)*sigma_y # generate errors for y (independent)
#'   M = beta_1[1] + X*beta_1[2]  + eps_M # generate latent mediator M
#'   y = beta_2[1] +  M*beta_2[2] + X*beta_2[3] + eps_Y # generate dependent variable
#'
#'   eps_m_star[,1]=rnorm(N)*sigma_m_star[1] # generate errors for m_star (independent)
#'   m_star[,1] =  M + eps_m_star[,1] # generate observed mediator indicators m_star
#'   if(m_ind>1){
#'     for(i in 2:(m_ind))   {
#'       eps_m_star[,i]=rnorm(N)*sigma_m_star[i] # generate errors for m_star (independent)
#'       m_star[,i] =  lambda[(i-1),1] + M*lambda[(i-1),2] + X*beta_2[-c(1,2)] + eps_m_star[,i]
#'     }
#'   }
#'
#'   eps_y_star[,1]=rnorm(N)*sigma_y_star[1] # generate errors for y_star (independent)
#'   y_star[,1] =  y + eps_y_star[,1] # generate observed dependent variable indicators y_star
#'   if(y_ind>1){
#'     for(i in 2:(y_ind)){
#'       eps_y_star[,i]=rnorm(N)*sigma_y_star[i] # generate errors for y_star (independent)
#'       y_star[,i] =  tau[(i-1),1] + y*tau[(i-1),2] + eps_y_star[,i]
#'     }
#'   }
#'   list(X = X, M = M, m_star = m_star, y = y, y_star=y_star)
#' }
#' m_ind = 2; y_ind = 2;
#' sigma_M = 1^.5 # error std M
#' sigma_y = 1^.5 # error std y
#' sigma_m_star = c(.3,.5)^.5 #c(1,2)^.5
#' sigma_y_star = c(.5,.3)^.5  #c(2,1)^.5
#' beta_1 = c(.5,1)
#' beta_2 = c(1,3,0)
#' lambda = matrix(c(1,1.5),ncol=2)
#' tau = matrix(c(1,2),ncol = 2)
#' k=length(beta_1)-1
#' nobs = 1000   # number of observations
#' X = runif(nobs) # generate random X from a uniform distribution
#' Data = SimMeasurementCont( beta_1, beta_2 , lambda, tau, m_ind, y_ind, sigma_M, sigma_m_star,
#'                           sigma_y, sigma_y_star, nobs, X)
#' R = 5000; burnin = 3000
#'
#' A_M=rep(100,2);
#' A_Y=c(100,100,1)
#'
#' #Estimation
#' out = MeasurementCont(Data = Data, Prior = list(A_M = A_M, A_Y = A_Y),R=5000, burnin = 3000)
#'
#' #Results
#' colMeans(out$beta_M)
#' colMeans(out$beta_Y)
#'
#' BFMeasurementCont = exp(BFSD(post = out , prior = A_Y[3],burnin = 0))
### Description MeasurementCont estimates a partial mediation model with multiple categorical indicators for the mediator
# and the dependent variable using Hamiltonian Monte Carlo (HMC) with Stan

### Arguments:
# Data  list(X, m_star, y_star)
# Prior list(A_M,A_Y)
# R
# burnin number of MCMC draws before the posterior is converged

### Details:
## Model:
# M = beta_0M + Xbeta_1 + U_M   (eq.1)
# Y = beta_0Y + Mbeta_2 + Xbeta_3 + U_Y  (eq.2)
# indicator equations:
# m*_1 = M + U_m*_1
# m*_2 = lambda_01 + M + U_m*_2
# ...
# m*_k = lambda_0k-1 + M + U_m*_k
# y*_1 = M + U_y*_1
# y*_2 = tau_01 + M + U_y*_2
# ...
# y*_l = tau_0l-1 + M + U_y*_l

## Data = list(X, m_star, y_star)
# X(N x 1) treatment variable vector
# m_star(N x M_ind) mediator indicators' matrix
# y_star(N x Y_ind) dependent variable indicators' matrix
# Prior = list(A_M,A_Y) [optional]
# A_M vector of coefficients' prior variances of eq.1 (def: rep(100,2))
# A_Y vector of coefficients' prior variances of eq.2 (def: c(100,100,1))
# R number of MCMC iterations (def:10000)
# burnin number of MCMC iterations to be discarded from the draws
### Value:
# beta_M(R X 2)  matrix of eq.1 coefficients' draws
# beta_Y(R X 3)  matrix of eq.2 coefficients' draws
# lambda (M_ind X 2 X R) array of mediator indicator coefficients' draws.
# tau (Y_ind X 2 X R) array of dependent variable indicator coefficients' draws.
# Each slice is one draw, where rows represent the indicator equation and columns are the coefficients
# All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively
# ssq_m_star(R X M_ind) Matrix of mediator indicator equations' coefficients' error variance draws
# ssq_y_star(R X Y_ind) Matrix of dependent variable indicator equation's coefficients' error variance draws
# mu_draw(R X 1) vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
# var_draw(R X 1) vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
MeasurementCont = function(Data, Prior, R, burnin){

  ############################################
  ## Arash Laghaie 2019
  ############################################

  if(missing(Prior))
  { A_M = rep(10,2); A_Y = c(10,10,1);}
  else
  {
    if(is.null(Prior$A_M)) {A_M = rep(10,2)}
    else {A_M = sqrt(Prior$A_M)}
    if(is.null(Prior$A_Y)) {A_Y = c(10,10,1)}
    else {A_Y = sqrt(Prior$A_Y)}
  }

  X = Data$X
  # M = Data$M
  # Y = Data$Y
  N = length(Data$X);

  m_star = Data$m_star
  y_star = as.matrix(Data$y_star)
  m_ind =  dim(m_star)[2];
  y_ind =  dim(y_star)[2];


  stanfit = rstan::sampling(object = stanmodels$Measurement_Multi,
                            data = list(n=N,M_ind = m_ind, Y_ind = y_ind,A_M=A_M, A_Y=A_Y, X = X, m_star=t(m_star),y_star=t(y_star)),
                            pars = c("beta_1","beta_2","beta_3","ssq_M","ssq_Y","ssq_m_star","ssq_y_star","beta_0_Y","beta_0_M","lambda","tau","M","Y"),
                            chains = 1,iter = R, warmup = burnin)

  post = rstan::extract(stanfit)

  #Running a Gibbs sampler regression of latent Y on latent M and X to store the MCMC draws for Bayes factor computation
  out = RuniregGibbsMulti(Data = list(y=post$Y, M = post$M, X=as.matrix(X)), Prior = list(ssq=1,A=diag(1/A_Y^2)), Mcmc = list(R=R-burnin,sigmasq=post$ssq_Y))

  iota = matrix(rep(c(0,1),R-burnin), ncol=2, byrow=T)  #first measurement equations' fixed parameters
  post$lambda = abind::abind(iota,post$lambda,along=2)
  post$tau = abind::abind(iota,post$tau,along=2)
  return(list(beta_M = cbind(post$beta_0_M,post$beta_1), beta_Y = cbind(post$beta_0_Y,post$beta_2,post$beta_3),
              ssq_M=post$ssq_M, ssq_Y = post$ssq_Y, Mdraw = post$M, Ydraw = post$Y,
              ssq_m_star = post$ssq_m_star, ssq_y_star = post$ssq_y_star,
              lambda = post$lambda, tau = post$tau,
              mu_draw = out$mubeta[,3], var_draw = out$IR[3,3,]))
}
