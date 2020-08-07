#' Sampler for Partial Mediation Model with Multiple Categorical Indicator for the Mediator
#' @description
#' Estimates a partial mediation model with multiple categorical indicator for the mediator and observed dependent variable using a mixture of Metropolis-Hastings and Gibbs sampling
#'
#' @usage MeasurementMCat(Data,Prior,R)
#'
#' @param Data list(X, m_star, Y)
#' @param Prior list(A_M,A_Y)
#' @param R number of MCMC iterations, default = 10000
#'
#' @details
#' ## Model
#'
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
#'
#'
#' ## Argument Details
#'
#' ## \code{Data = list(X, m_star, Y)}
#' \describe{
#'   \item{X(N x 1)}{treatment variable vector}
#'   \item{m_star(N x M_ind)}{mediator indicators' matrix}
#'   \item{Y(N x 1)}{dependent variable vector}
#' }
#'
#'
#' ## \code{Prior = list(A_M,A_Y)} \[optional\]
#' \describe{
#'   \item{A_M}{vector of coefficients' prior variances of eq.1, default = rep(100,2)}
#'   \item{A_Y}{vector of coefficients' prior variances of eq.2, default = c(100,100,1)}
#' }
#'
#' @return
#' \describe{
#'   \item{beta_M(R X 2)}{matrix of eq.1 coefficients' draws}
#'   \item{beta_Y(R X 3)}{matrix of eq.2 coefficients' draws}
#'   \item{lambda (M_ind X 2 X R)}{array of indicator coefficients' draws. Each slice is one draw, where rows represent the indicator equation and columns are the coefficients. All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively.}
#'   \item{ssq_m_star(R X M_ind)}{Matrix of indicator equations' coefficients' error variance draws}
#'   \item{ssq_Y(R X 1)}{vector of eq.2 error variance draws}
#'   \item{mu_draw}{vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#'   \item{var_draw}{vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#' }
#'
#' @export
#' @examples
#' SimMeasurementMCat = function(X, beta_M, cutoff_M, beta_Y, Sigma_Y, M_ind, lambda, ssq_m_star){
#'
#'   nobs = dim(X)[1]
#'   m_star = m_tilde = matrix(double(nobs*M_ind), ncol = M_ind)
#'
#'   M = beta_M[1] + beta_M[2] * X + rnorm(nobs)  #cbind(rep(1,nobs),X)%*%beta_M + rnorm(nobs)
#'
#'   for(i in 1: M_ind){
#'     m_star[,i] = lambda[i] + M + sqrt(ssq_m_star[i])*rnorm(nobs);
#'     m_tilde[,i] = cut(m_star[,i], br = cutoff_M[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
#'   }
#'
#'   Y = beta_Y[1] + beta_Y[2] * M + beta_Y[3] * X + rnorm(nobs)*Sigma_Y
#'                                                      #cbind(rep(1,nobs),cbind(M,X))%*%beta_Y + rnorm(nobs)
#'
#'   return(list(Y = Y, M = M, m_tilde = m_tilde, X = X,
#'               beta_M = beta_M, beta_Y = beta_Y,
#'               lambda = lambda, ssq_m_star = ssq_m_star, m_star = m_star, cutoff_M = cutoff_M,
#'               k_M=dim(cutoff_M)[2]-1, M_ind=M_ind))
#' }
#'
#' M_ind = 2
#' Mcut =  8
#' nobs= 1000
#' X=as.matrix(runif(nobs,min=0, max=1))
#' beta_M = c(.5,1)
#' beta_Y = c(1, 2, 0)
#' Sigma_Y = 1^.5
#' ssq_m_star = c(.5,.7)
#' lambda = c(0,-.5)   #the intercepts for the latent M indicators w. measurement
#'                     #error (first intercept should always be 0)
#'
#' cutoff_M = matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
#'                     -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Mcut, byrow = T)
#' DataMCat = SimMeasurementMCat(X, beta_M, cutoff_M, beta_Y, Sigma_Y, M_ind, lambda, ssq_m_star)
#'
#' #estimation
#' Data = list(X=cbind(rep(1,length(DataMCat$X)),DataMCat$X), m_tilde=as.matrix(DataMCat$m_tilde),
#'             Y= as.matrix(DataMCat$Y) ,k=Mcut-1, M_ind=dim(DataMCat$m_tilde)[2])
#' out = MeasurementMCat(Data=Data, R=R)
#'
#'
#' #results
#' colMeans(out$beta_M)
#' colMeans(out$beta_Y)
#' apply(out$cutoff_M,c(1,2),FUN = mean)
### Description MeasurementMCat estimates a partial mediation model with multiple categorical indicator for the mediator
# and observed dependent variable using a mixture of Metropolis-Hastings and Gibbs sampling
### Arguments:
# Data  list(X, m_tilde, Y)
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
## Data = list(X, m_tilde, Y)
# X(N x 1) treatment variable vector
# m_tilde(N x M_ind) mediator indicators' matrix
# Y(N x 1) dependent variable vector
# Prior = list(A_M,A_Y) [optional]
# A_M vector of coefficients' prior variances of eq.1 (def: rep(100,2))
# A_Y vector of coefficients' prior variances of eq.2 (def: c(100,100,1))
# R number of MCMC iterations (def:10000)
### Value:
# beta_M(R X 2)  matrix of eq.1 coefficients' draws
# beta_Y(R X 3)  matrix of eq.2 coefficients' draws
# lambda (M_ind X 2 X R) array of indicator coefficients' draws.
# Each slice is one draw, where rows represent the indicator equation and columns are the coefficients
# All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively
# ssq_m_star(R X M_ind) Matrix of indicator equations' coefficients' error variance draws
# ssq_Y(R X 1) vector of eq.2 error variance draws
# mu_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
# var_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
MeasurementMCat=function(Data,Prior,R=10000){
  ############################################
  # Arash Laghaie 2019
  # The ordered probit part of the inference and the R shell is based on rordprobitGibbs {bayesm}
  ############################################
  # define functions needed
  #  dstartoc is a fuction to transfer dstar to its cut-off value
  dstartoc=function(dstar) {c(-100, 0, cumsum(exp(dstar)), 100)}

  # compute conditional likelihood of data given cut-offs
  #
  lldstar=function(dstar,y,mu){
    gamma=dstartoc(dstar)
    arg = pnorm(gamma[y+1]-mu)-pnorm(gamma[y]-mu)
    epsilon=1.0e-50
    arg=ifelse(arg < epsilon,epsilon,arg)
    return(sum(log(arg)))
  }
  #
  # ----------------------------------------------------------------------
  #
  # check arguments
  #
  if(missing(Data)) {pandterm("Requires Data argument -- list of y and X")}
  if(is.null(Data$X)) {pandterm("Requires Data element X")}
  X=Data$X
  if(is.null(Data$m_tilde)) {pandterm("Requires Data element m_tilde")}
  y=Data$m_tilde
  if(is.null(Data$Y)) {pandterm("Requires Data element Y")}
  dep=Data$Y
  # if(is.null(Data$beta_2)) {pandterm("Requires Data element beta_2")}
  # beta_2=Data$beta_2
  if(is.null(Data$k)) {pandterm("Requires Data element k")}
  k=Data$k
  if(is.null(Data$M_ind)) {pandterm("Requires Data element M_ind")}
  M_ind=Data$M_ind


  nvar=ncol(X)
  nobs=dim(y)[1]
  ndstar = k-2         # number of dstar being estimated
  ncuts = k+1          # number of cut-offs (including zero and two ends)
  ncut = ncuts-3       # number of cut-offs being estimated c[1]=-100, c[2]=0, c[k+1]=100


  # cutoff_Y_init = Data$cutoff_M_init
  # Y_tilde_init = Data$M_tilde_init
  # ssq_y_tilde_init = Data$ssq_y_tilde_init
  # beta_init = Data$beta_init
  # beta_2_init = Data$beta_2_init
  # Y_init = Data$M_init
  # beta_tilde_init = Data$beta_tilde_init

  #
  # check data for validity
  #
  if(dim(y)[1] != nrow(X) ) {pandterm("y and X not of same row dim")}
  if(  sum(unique(y[,1]) %in% (1:k) ) < length(unique(y[,1])) )     #Tests only y_star1 but I really don't care! :D
  {pandterm("some value of y is not vaild")}

  #
  # check for Prior
  #
  if(missing(Prior))
  { betabar=c(rep(0,nvar)); A=diag(1/rep(100,2)); Ad=diag(ndstar); dstarbar=c(rep(0,ndstar));
  betabar_2=c(rep(0,nvar+1)); A_2=diag(1/c(100,100,1));}
  else
  {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))}
    else {betabar=Prior$betabar}
    if(is.null(Prior$A_M)) {A=diag(1/rep(100,2))}
    else {A=diag(1/Prior$A_M)}                       #A (A_M) prior variance of beta_M vector
    if(is.null(Prior$betabar_2)) {betabar_2=c(rep(0,nvar+1))}
    else {betabar_2=Prior$betabar_2}
    if(is.null(Prior$A_Y)) {A_2=diag(1/c(100,100,1))}
    else {A_2=diag(1/Prior$A_Y)}                     #A_2 (A_Y) prior variance of beta_Y vector
    if(is.null(Prior$Ad)) {Ad=diag(ndstar)}
    else {Ad=Prior$Ad}
    if(is.null(Prior$dstarbar)) {dstarbar=c(rep(0,ndstar))}
    else {dstarbar=Prior$dstarbar}
  }
  #
  # check dimensions of Priors
  #

  if(ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar)
  {pandterm(paste("bad dimensions for A",dim(A)))}
  if(length(betabar) != nvar)
  {pandterm(paste("betabar wrong length, length= ",length(betabar)))}
  if(ncol(A_2) != nrow(A_2) || ncol(A_2) != nvar+1 || nrow(A_2) != nvar+1)
  {pandterm(paste("bad dimensions for A_2",dim(A_2)))}
  if(length(betabar_2) != nvar+1)
  {pandterm(paste("betabar_2 wrong length, length= ",length(betabar_2)))}
  if(ncol(Ad) != nrow(Ad) || ncol(Ad) != ndstar || nrow(Ad) != ndstar)
  {pandterm(paste("bad dimensions for Ad",dim(Ad)))}
  if(length(dstarbar) != ndstar)
  {pandterm(paste("dstarbar wrong length, length= ",length(dstarbar)))}

  #
  # check MCMC argument
  #
  keep = 1
  nprint = 100
  s = 2.38/sqrt(ndstar)
  # if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
  # else
  # {
  #   if(is.null(Mcmc$R))
  #   {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
  #   if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
  #   if(is.null(Mcmc$nprint)) {nprint=100} else {nprint=Mcmc$nprint}
  #   if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
  #   if(is.null(Mcmc$s)) {s=2.38/sqrt(ndstar)} else {s=Mcmc$s} #2.38 is the RRscaling
  # }
  #
  # print out problem
  #
  cat(" ", fill=TRUE)
  cat("Starting Gibbs Sampler for mediation LVM Model",fill=TRUE)
  cat("   with ",nobs,"observations",fill=TRUE)
  cat(" ", fill=TRUE)
  cat("Table of m_tilde values",fill=TRUE)
  for(i in 1:M_ind) print(table(y[,i]))
  cat(" ",fill=TRUE)
  # cat("Prior Parms: ",fill=TRUE)
  # cat("betabar",fill=TRUE)
  # print(betabar)
  # cat(" ", fill=TRUE)
  # cat("A",fill=TRUE)
  # print(A)
  # cat(" ", fill=TRUE)
  # cat("betabar_2",fill=TRUE)
  # print(betabar_2)
  # cat(" ", fill=TRUE)
  # cat("A_2",fill=TRUE)
  # print(A_2)
  # cat(" ", fill=TRUE)
  # cat("dstarbar",fill=TRUE)
  # print(dstarbar)
  # cat(" ", fill=TRUE)
  # cat("Ad",fill=TRUE)
  # print(Ad)
  # cat(" ", fill=TRUE)
  # cat("MCMC parms: ",fill=TRUE)
  # cat("R= ",R," keep= ",keep," nprint= ",nprint,"s= ",s, fill=TRUE)
  # cat(" ",fill=TRUE)

  # use (-Hessian+Ad)^(-1) evaluated at betahat as the basis of the
  # covariance matrix for the random walk Metropolis increments

  betahat = chol2inv(chol(crossprod(X,X)))%*% crossprod(X,rowMeans(y))    ###betahat is computed based on the average of y_stars
  dstarini = c(cumsum(c( rep(0.1, ndstar))))     # set initial value for dstar
  dstarout = optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=X%*%betahat, y=rowMeans(y))
  inc.root=chol(chol2inv(chol((-dstarout$hessian+Ad))))  # chol((H+Ad)^-1)

  ###################################################################
  draws=MeasurementMCatCpp(dep,y,X,k,A,betabar,Ad,A_2,betabar_2,
                                            s,inc.root,dstarbar,betahat,
                                            M_ind,
                                            R,keep,nprint)
                                            # cutoff_Y_init, Y_tilde_init, beta_tilde_init, ssq_y_tilde_init, beta_init, beta_2_init, Y_init)

  ###################################################################

  draws$cutoff_M=draws$cutoff_M[,2:k,]


  return(draws)
}
