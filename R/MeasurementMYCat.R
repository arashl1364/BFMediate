#' Sampler for Partial Mediation Model with Multiple Categorical Indicator for the Mediator and DV
#' @description
#' Estimates a partial mediation model with multiple categorical indicator for the mediator and the dependent variable using a mixture of Metropolis-Hastings and Gibbs sampling
#'
#' @usage MeasurementMYCat(Data,Prior,R)
#'
#' @param Data list(X, m_tilde, y_tilde)
#' @param Prior list(A_M,A_Y)
#' @param R number of MCMC iterations, default = 10000
#'
#' @details
#' ## Model:
#' (eq.1) \deqn{M = \beta_{0M} + X\beta_1 + U_M}{M = \beta_0M + X\beta_1 + U_M}
#' (eq.2) \deqn{Y = \beta_{0Y} + M\beta_2 + X\beta_3 + U_Y}{Y = \beta_0Y + M\beta_2 + X\beta_3 + U_Y}
#'
#' ## Indicator equations:
#' \deqn{m^*_1 = M + U_{m^*_1}}{m*_1 = M + U_{m*_1}}
#' \deqn{\tilde{m}_1 = OrdProbit(m^*_1 ,C_{m_1})}{˜m_1 = OrdProbit(m*_1,C_{m_1})}
#' \deqn{m^*_2 = \lambda_{01} + M + U_{m^*_2}}{m*_2 = \lambda_01 + M + U_{m*_2}}
#' \deqn{\tilde{m}_2 = OrdProbit(m^*_2 ,C_{m_2})}{˜m_2 = OrdProbit(m*_2,C_{m_2})}
#' \deqn{...}
#' \deqn{m^*_k = \lambda_{0k-1} + M + U_{m^*_k}}{m*_k = \lambda_0k-1 + M + U_{m*_k}}
#' \deqn{\tilde{m}_k = OrdProbit(m^*_k ,C_{m_k})}{˜m_k = OrdProbit(m*_k,C_{m_k})}
#' \deqn{y^*_1 = M + U_{y^*_1}}{y*_1 = M + U_{y*_1}}
#' \deqn{\tilde{y}_1 = OrdProbit(y^*_1 ,C_{y_1})}{˜y_1 = OrdProbit(y*_1,C_{y_1})}
#' \deqn{y^*_2 = \tau_{01} + M + U_{y^*_2}}{y*_2 = \tau_01 + M + U_{y*_2}}
#' \deqn{\tilde{y}_2 = OrdProbit(y^*_2 ,C_{y_2})}{˜y_2 = OrdProbit(y*_2,C_{y_2})}
#' \deqn{...}
#' \deqn{y^*_l = \tau_{0l-1} + M + U_{y^*_l}}{y*_l = \tau_0l-1 + M + U_{y*_l}}
#' \deqn{\tilde{y}_l = OrdProbit(y^*_l ,C_{y_l})}{˜y_l = OrdProbit(y*_l,C_{y_l})}
#'
#'
#' ## Argument Details
#'
#' ## \code{Data = list(X, m_tilde, y_tilde)}
#' \describe{
#'   \item{X(N x 1)}{treatment variable vector}
#'   \item{m_tilde(N x M_ind)}{mediator indicators' matrix}
#'   \item{y_tilde(N x Y_ind)}{dependent variable indicators' matrix}
#' }
#'
#' ## \code{Prior = list(A_M,A_Y)} \[optional\]
#' \describe{
#'   \item{A_M}{vector of coefficients' prior variances of eq.1, default = rep(100,2)}
#'   \item{A_Y}{vector of coefficients' prior variances of eq.2, default = c(100,100,1)}
#' }
#'
#' @return
#' \describe{
#'   \item{beta_1(R X 2)}{matrix of eq.1 coefficients' draws}
#'   \item{beta_2(R X 3)}{matrix of eq.2 coefficients' draws}
#'   \item{lambda (M_ind X 2 X R)}{array of mediator indicator coefficients' draws. Each slice is one draw, where rows represent the indicator equation and columns are the coefficients. All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively.}
#'   \item{ssq_m_star(R X M_ind)}{Matrix of mediator indicator equations' coefficients' error variance draws}
#'   \item{ssq_y_star(R X Y_ind)}{Matrix of dependent variable indicator equations' coefficients' error variance draws}
#'   \item{mu_draw}{vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#'   \item{var_draw}{vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#' }
#'
#' @export
#' @examples
#' SimMeasurementMYCat = function(X, beta_1, cutoff_M, beta_2, cutoff_Y, M_ind, Y_ind,
#'                                lambda, tau, ssq_m_star, ssq_y_star){
#'
#'   nobs = dim(X)[1]
#'   m_tilde = m_star = matrix(double(nobs*M_ind), ncol = M_ind)
#'   y_tilde = y_star = matrix(double(nobs*Y_ind), ncol = Y_ind)
#'
#'   M = cbind(rep(1,nobs),X)%*%beta_1 + rnorm(nobs)
#'   for(i in 1: M_ind){
#'     m_star[,i] = lambda[i] + M + sqrt(ssq_m_star[i])*rnorm(nobs);
#'     m_tilde[,i] = cut(m_star[,i], br = cutoff_M[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
#'   }
#'
#'   Y = cbind(rep(1,nobs),cbind(M,X))%*%beta_2 + rnorm(nobs)
#'   for(i in 1: Y_ind){
#'     y_star[,i] = tau[i] + Y + sqrt(ssq_y_star[i])*rnorm(nobs);
#'     y_tilde[,i] = cut(y_star[,i], br = cutoff_Y[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
#'   }
#'
#'   return(list(Y = Y, M = M, y_tilde = y_tilde, m_tilde = m_tilde, X = X,
#'               k_M=dim(cutoff_M)[2]-1, beta_1 = beta_1, lambda = lambda,
#'               ssq_m_star = ssq_m_star, m_star = m_star, cutoff_M = cutoff_M,
#'               k_Y=dim(cutoff_Y)[2]-1, beta_2 = beta_2, tau = tau, ssq_y_star = ssq_y_star,
#'               y_star = y_star, cutoff_Y = cutoff_Y, M_ind=M_ind, Y_ind=Y_ind))
#' }
#'
#' M_ind = 2
#' Y_ind = 2
#' Mcut = Ycut = 8
#' nobs=2000
#' X=as.matrix(runif(nobs,min=0, max=1))
#' beta_1 = c(.5,1)
#' beta_2 = c(.7, 1.5, 0)
#' ssq_m_star = c(.5,.3)
#' lambda = c(0,-.5) #the intercepts for the latent M indicators w. measurement error
#' ssq_y_star = c(.2,.2)
#' tau = c(0,-.5)
#' cutoff_M = matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
#'                     -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Mcut, byrow = T)
#' cutoff_Y =  matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
#'                      -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Ycut, byrow = T)
#' DataMYCat = SimMeasurementMYCat(X, beta_1, cutoff_M, beta_2, cutoff_Y, M_ind,
#'                                 Y_ind, lambda, tau, ssq_m_star, ssq_y_star)
#'
#'
#' #estimation
#' Mcut = max(DataMYCat$m_tilde) + 1
#' Ycut = max(DataMYCat$y_tilde) + 1
#' Data = list(X=cbind(rep(1,length(DataMYCat$X)),DataMYCat$X), m_tilde=as.matrix(DataMYCat$m_tilde),
#'             y_tilde=as.matrix(DataMYCat$y_tilde), k_M = Mcut-1, k_Y=Ycut-1,
#'             M_ind=dim(as.matrix(DataMYCat$m_tilde))[2], Y_ind=dim(as.matrix(DataMYCat$y_tilde))[2])
#' out = MeasurementMYCat(Data=Data,R = 10000)
#'
#' #results
#' colMeans(out$beta_1)
#' colMeans(out$beta_2)
#' apply(out$lambdadraw,MARGIN = c(1,2),FUN = mean)
#' apply(out$taudraw,MARGIN = c(1,2),FUN = mean)
#'
#' apply(out$cutoff_M,c(1,2),FUN = mean)
#' apply(out$cutoff_Y,c(1,2),FUN = mean)
### Description MeasurementMYCat estimates a partial mediation model with multiple categorical indicator for the mediator
# and the dependent variable using a mixture of Metropolis-Hastings and Gibbs sampling

### Arguments:
# Data  list(X, m_tilde, y_tilde)
# Prior list(A_M,A_Y)
# R

### Details:
## Model:
# M = beta_0M + Xbeta_1 + U_M   (eq.1)
# Y = beta_0Y + Mbeta_2 + Xbeta_3 + U_Y  (eq.2)
# indicator equations:
# m*_1 = M + U_m*_1
# ˜m_1 = OrdProbit(m*_1,C_m_1)
# m*_2 = lambda_01 + M + U_m*_2
# ˜m_2 = OrdProbit(m*_2,C_m_2)
# ...
# m*_k = lambda_0k-1 + M + U_m*_k
# ˜m_k = OrdProbit(m*_k,C_m_k)
# y*_1 = M + U_y*_1
# ˜y_1 = OrdProbit(y*_1,C_y_1)
# y*_2 = tau_01 + M + U_y*_2
# ˜y_2 = OrdProbit(y*_2,C_y_2)
# ...
# y*_l = tau_0l-1 + M + U_y*_l
# ˜y_l = OrdProbit(y*_l,C_y_l)

## Data = list(X, m_tilde, y_tilde)
# X(N x 1) treatment variable vector
# m_tilde(N x M_ind) mediator indicators' matrix
# y_tilde(N x Y_ind) dependent variable indicators' matrix
# Prior = list(A_M,A_Y) [optional]
# A_M vector of coefficients' prior variances of eq.1 (def: rep(100,2))
# A_Y vector of coefficients' prior variances of eq.2 (def: c(100,100,1))
# R number of MCMC iterations (def:10000)

### Value:
# beta_1(R X 2)  matrix of eq.1 coefficients' draws
# beta_2(R X 3)  matrix of eq.2 coefficients' draws
# lambda (M_ind X 2 X R) array of mediator indicator coefficients' draws.
# tau (Y_ind X 2 X R) array of dependent variable indicator coefficients' draws.
# Each slice is one draw, where rows represent the indicator equation and columns are the coefficients
# All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively
# ssq_m_star(R X M_ind) Matrix of mediator indicator equations' coefficients' error variance draws
# ssq_y_star(R X Y_ind) Matrix of dependent variable indicator equations' coefficients' error variance draws
# mu_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
# var_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)

MeasurementMYCat=function(Data,Prior,R=10000){    #,Mcmc){
  # Rcpp::sourceCpp('Mediation_Ordered_Multi_Merr.cpp')
  #
  # revision history:
  #   3/07  Hsiu-Wen Liu
  #   3/07  fixed naming of dstardraw rossi
  #   19/10/2018 Arash Laghaie
  # purpose:
  #   draw from posterior for ordered probit using Gibbs Sampler
  #   and metropolis RW
  #
  # Arguments:
  #   Data - list of X,y_tilde,k_Y
  #     X is nobs x nvar_Y, y_star is nobs vector of 1,2,.,k_Y (ordinal variable)
  #   Prior - list of A_Y, beta_2_bar
  #     A_Y is nvar_Y x nvar_Y prior preci matrix
  #     beta_2_bar is nvar_Y x 1 prior mean
  #     Ad_Y is ndstar_Y x ndstar_Y prior preci matrix of dstar_Y (ncut_Y is number of cut-offs being estimated)
  #     dstarbar_Y is ndstar_Y x 1 prior mean of dstar_Y
  #   Mcmc
  #     R is number of draws
  #     keep is thinning parameter
  #     nprint - print estimated time remaining on every nprint'th draw
  #     s_Y is scale parameter of random work Metropolis
  #
  # Output:
  #   list of beta_2_draws and cutdraws
  #
  # Model:
  #    z=Xbeta_2 + e  < 0  e ~N(0,1)
  #    y_star=1,..,k_Y, if z~c(c[k_Y], c[k_Y+1])
  #
  #    cutoffs = c[1],..,c[k_Y+1]
  #    dstar_Y = dstar_Y[1],dstar_Y[k_Y-2]
  #    set c[1]=-100, c[2]=0, ...,c[k_Y+1]=100
  #
  #    c[3]=exp(dstar_Y[1]),c[4]=c[3]+exp(dstar_Y[2]),...,
  #    c[k_Y]=c[k_Y-1]+exp(datsr[k_Y-2])
  #
  # Note: 1. length of dstar_Y = length of cutoffs - 3
  #       2. Be careful in assessing prior parameter, Ad_Y.  .1 is too small for many applications.
  #
  # Prior: beta_2 ~ N(beta_2_bar,A_Y^-1)
  #        dstar_Y ~ N(dstarbar_Y, Ad_Y^-1)
  #
  #
  # ----------------------------------------------------------------------
  # Rcpp::sourceCpp('Mediation_Ordered_Multi_Merr.cpp')
  # define functions needed
  #  dstartoc is a fuction to transfer dstar_Y to its cut-off value

  dstartoc=function(dstar) {c(-100, 0, cumsum(exp(dstar)), 100)}

  # compute conditional likelihood of data given cut-offs
  #
  lldstar=function(dstar_Y,y_tilde,mu){
    gamma=dstartoc(dstar_Y)
    arg = pnorm(gamma[y_tilde+1]-mu)-pnorm(gamma[y_tilde]-mu)
    epsilon=1.0e-50
    arg=ifelse(arg < epsilon,epsilon,arg)
    return(sum(log(arg)))
  }
  #
  # ----------------------------------------------------------------------
  #
  # check arguments
  #
  if(missing(Data)) {pandterm("Requires Data argument -- list of y_tilde and X")}
  if(is.null(Data$X)) {pandterm("Requires Data element X")}
  X=Data$X
  if(is.null(Data$y_tilde)) {pandterm("Requires Data element y_tilde")}
  y_tilde=Data$y_tilde
  if(is.null(Data$k_Y)) {pandterm("Requires Data element k_Y")}
  k_Y=Data$k_Y
  if(is.null(Data$m_tilde)) {pandterm("Requires Data element m_tilde")}
  m_tilde=Data$m_tilde
  if(is.null(Data$k_M)) {pandterm("Requires Data element k_M")}
  k_M=Data$k_M
  if(is.null(Data$M_ind)) {pandterm("Requires Data element M_ind")}
  M_ind=Data$M_ind
  if(is.null(Data$Y_ind)) {pandterm("Requires Data element Y_ind")}
  Y_ind=Data$Y_ind

  # #fixing parameters for testing the sampler
  # cutoffs_M_init = Data$cutoffs_M_init
  # M_tilde_init = Data$M_tilde_init
  # beta_m_tilde_init = Data$beta_m_tilde_init
  # ssq_m_tilde_init = Data$ssq_m_tilde_init
  # beta_init = Data$beta_init
  # M_init = Data$M_init
  #
  # cutoffs_Y_init = Data$cutoffs_Y_init
  # Y_tilde_init = Data$Y_tilde_init
  # beta_y_tilde_init = Data$beta_y_tilde_init
  # ssq_y_tilde_init = Data$ssq_y_tilde_init
  # beta_2_init = Data$beta_2_init
  # Y_init = Data$Y_init

  # cutoffs_M_init,
  # M_tilde_init, beta_tilde_init, ssq_m_tilde_init, beta_init, M_init,
  # cutoffs_Y_init,
  # Y_tilde_init, beta_2_tilde_init, ssq_y_tilde_init, beta_2_init, Y_init)



  nvar_Y=ncol(X)+1
  nobs=dim(y_tilde)[1]
  ndstar_Y = k_Y-2         # number of dstar_Y being estimated
  ncuts_Y = k_Y+1          # number of cut-offs (including zero and two ends)
  ncut_Y = ncuts_Y-3       # number of cut-offs being estimated c[1]=-100, c[2]=0, c[k_Y+1]=100
  nvar_M=ncol(X)
  ndstar_M = k_M-2         # number of dstar_M being estimated
  ncuts_M = k_M+1          # number of cut-offs (including zero and two ends)
  ncut_M = ncuts_M-3       # number of cut-offs being estimated c[1]=-100, c[2]=0, c[k_M+1]=100


  #
  # check data for validity
  #
  if(dim(y_tilde)[1] != nrow(X) ) {pandterm("y_tilde and X not of same row dim")}
  if(  sum(unique(y_tilde[,1]) %in% (1:k_Y) ) < length(unique(y_tilde[,1])) )     #Tests only y_star1 but I really don't care! :D
  {pandterm("some value of y_tilde is not vaild")}
  if(dim(m_tilde)[1] != nrow(X) ) {pandterm("m_tilde and X not of same row dim")}
  if(  sum(unique(m_tilde[,1]) %in% (1:k_M) ) < length(unique(m_tilde[,1])) )     #Tests only m_star1 but I really don't care! :D
  {pandterm("some value of m_tilde is not vaild")}

  #
  # check for Prior
  #
  if(missing(Prior))
  { betabar=c(rep(0,nvar_M)); A_M=diag(1/rep(100,2)); Ad_M=diag(ndstar_M); dstarbar_M=c(rep(0,ndstar_M))
  beta_2_bar=c(rep(0,nvar_Y)); A_Y=diag(1/c(100,100,1)); Ad_Y=diag(ndstar_Y); dstarbar_Y=c(rep(0,ndstar_Y))}
  else
  {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar_M))}
    else {betabar=Prior$betabar}
    if(is.null(Prior$A_M)) {A_M=diag(1/rep(100,2))}
    else {A_M=diag(1/Prior$A_M)}                     #A_M prior variance of beta_1 vector
    if(is.null(Prior$Ad_M)) {Ad_M=diag(ndstar_M)}
    else {Ad_M=Prior$Ad_M}
    if(is.null(Prior$dstarbar_M)) {dstarbar_M=c(rep(0,ndstar_M))}
    else {dstarbar_M=Prior$dstarbar_M}

    if(is.null(Prior$beta_2_bar)) {beta_2_bar=c(rep(0,nvar_Y))}
    else {beta_2_bar=Prior$beta_2_bar}
    if(is.null(Prior$A_Y)) {A_Y=diag(1/c(100,100,1))}   #A_Y prior variance of beta_2 vector
    else {A_Y=diag(1/Prior$A_Y)}
    if(is.null(Prior$Ad_Y)) {Ad_Y=diag(ndstar_Y)}
    else {Ad_Y=Prior$Ad_Y}
    if(is.null(Prior$dstarbar_Y)) {dstarbar_Y=c(rep(0,ndstar_Y))}
    else {dstarbar_Y=Prior$dstarbar_Y}
  }
  #
  # check dimensions of Priors
  #

  if(ncol(A_M) != nrow(A_M) || ncol(A_M) != nvar_M || nrow(A_M) != nvar_M)
  {pandterm(paste("bad dimensions for A_M",dim(A_M)))}
  if(length(betabar) != nvar_M)
  {pandterm(paste("betabar wrong length, length= ",length(betabar)))}
  if(ncol(Ad_M) != nrow(Ad_M) || ncol(Ad_M) != ndstar_M || nrow(Ad_M) != ndstar_M)
  {pandterm(paste("bad dimensions for Ad_M",dim(Ad_M)))}
  if(length(dstarbar_M) != ndstar_M)
  {pandterm(paste("dstarbar_M wrong length, length= ",length(dstarbar_M)))}

  if(ncol(A_Y) != nrow(A_Y) || ncol(A_Y) != nvar_Y || nrow(A_Y) != nvar_Y)
  {pandterm(paste("bad dimensions for A_Y",dim(A_Y)))}
  if(length(beta_2_bar) != nvar_Y)
  {pandterm(paste("beta_2_bar wrong length, length= ",length(beta_2_bar)))}
  if(ncol(Ad_Y) != nrow(Ad_Y) || ncol(Ad_Y) != ndstar_Y || nrow(Ad_Y) != ndstar_Y)
  {pandterm(paste("bad dimensions for Ad_Y",dim(Ad_Y)))}
  if(length(dstarbar_Y) != ndstar_Y)
  {pandterm(paste("dstarbar_Y wrong length, length= ",length(dstarbar_Y)))}

  #
  # check MCMC argument
  #
  keep = 1
  nprint = 100
  s_M = 2.38/sqrt(ndstar_M)
  s_Y = 2.38/sqrt(ndstar_Y)
  # if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
  # else
  # {
  #   if(is.null(Mcmc$R))
  #   {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
  #   if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
  #   if(is.null(Mcmc$nprint)) {nprint=100} else {nprint=Mcmc$nprint}
  #   if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
  #   if(is.null(Mcmc$s_M)) {s_M=2.38/sqrt(ndstar_M)} else {s_M=Mcmc$s_M} #2.38 is RRscaling
  #   if(is.null(Mcmc$s_Y)) {s_Y=2.38/sqrt(ndstar_Y)} else {s_Y=Mcmc$s_Y}
  # }
  #
  # print out problem
  #
  cat(" ", fill=TRUE)
  cat("Starting Gibbs Sampler for Ordered Probit Model",fill=TRUE)
  cat("   with ",nobs,"observations",fill=TRUE)
  cat(" ", fill=TRUE)
  cat("Table of y_tilde values",fill=TRUE)
  for(i in 1:Y_ind) print(table(y_tilde[,i]))
  cat("Table of m_tilde values",fill=TRUE)
  for(i in 1:M_ind) print(table(m_tilde[,i]))
  cat(" ",fill=TRUE)
  cat("Prior Parms: ",fill=TRUE)
  cat("betabar",fill=TRUE)
  print(betabar)
  cat("beta_2_bar",fill=TRUE)
  print(beta_2_bar)
  cat(" ", fill=TRUE)
  cat("A_M",fill=TRUE)
  print(A_M)
  cat("A_Y",fill=TRUE)
  print(A_Y)
  cat(" ", fill=TRUE)
  cat("dstarbar_M",fill=TRUE)
  print(dstarbar_M)
  cat("dstarbar_Y",fill=TRUE)
  print(dstarbar_Y)
  cat(" ", fill=TRUE)
  cat("Ad_M",fill=TRUE)
  print(Ad_M)
  cat("Ad_Y",fill=TRUE)
  print(Ad_Y)
  cat(" ", fill=TRUE)
  cat("MCMC parms: ",fill=TRUE)
  cat("R= ",R," keep= ",keep," nprint= ",nprint,"s_M= ",s_M,"s_Y= ",s_Y, fill=TRUE)
  cat(" ",fill=TRUE)

  # use (-Hessian+Ad_M)^(-1) evaluated at betahat as the basis of the
  # covariance matrix for the random walk Metropolis increments

  betahat = chol2inv(chol(crossprod(X,X)))%*% crossprod(X,rowMeans(m_tilde))
  dstarini = c(cumsum(c( rep(0.1, ndstar_M))))     # set initial value for dstar_M
  dstarout = optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=X%*%betahat, y=rowMeans(m_tilde))
  inc.root_M=chol(chol2inv(chol((-dstarout$hessian+Ad_M))))  # chol((H+Ad_Y)^-1)


  # use (-Hessian+Ad_Y)^(-1) evaluated at beta_2_hat as the basis of the
  # covariance matrix for the random walk Metropolis increments

  Xmstar=cbind(X[,1],rowMeans(m_tilde),X[,2])
  beta_2_hat = chol2inv(chol(crossprod(Xmstar,Xmstar)))%*% crossprod(Xmstar,rowMeans(y_tilde))
  dstarini = c(cumsum(c( rep(0.1, ndstar_Y))))     # set initial value for dstar_Y
  dstarout = optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=Xmstar%*%beta_2_hat, y=rowMeans(y_tilde))
  inc.root_Y=chol(chol2inv(chol((-dstarout$hessian+Ad_Y))))  # chol((H+Ad_Y)^-1)

  ###################################################################
  ###################################################################
  # Keunwoo Kim
  # 08/20/2014
  # Modified by Arash Laghaie
  # 12/27/2018
  ###################################################################
  draws= MeasurementMYCatCpp(X, m_tilde, y_tilde, k_M, k_Y,M_ind,Y_ind,
                                     A_M, betabar, Ad_M, s_M, inc.root_M, dstarbar_M, betahat,
                                     A_Y, beta_2_bar, Ad_Y, s_Y, inc.root_Y, dstarbar_Y, beta_2_hat,
                                     R, keep, nprint)
  ###################################################################

  draws$cutoff_M=draws$cutoff_M[,2:k_M,]
  # attributes(draws$cutdraw_M)$class="bayesm.mat"
  # attributes(draws$betadraw)$class="bayesm.mat"
  # attributes(draws$dstardraw_M)$class="bayesm.mat"
  # attributes(draws$cutdraw_M)$mcpar=c(1,R,keep)
  # attributes(draws$betadraw)$mcpar=c(1,R,keep)
  # attributes(draws$dstardraw_M)$mcpar=c(1,R,keep)

  draws$cutoff_Y=draws$cutoff_Y[,2:k_Y,]
  # attributes(draws$cutdraw_Y)$class="bayesm.mat"
  # attributes(draws$beta_2_draw)$class="bayesm.mat"
  # attributes(draws$dstardraw_Y)$class="bayesm.mat"
  # attributes(draws$cutdraw_Y)$mcpar=c(1,R,keep)
  # attributes(draws$beta_2_draw)$mcpar=c(1,R,keep)
  # attributes(draws$dstardraw_Y)$mcpar=c(1,R,keep)

  return(draws)
}
