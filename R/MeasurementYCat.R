#' Estimates a partial mediation model with multiple categorical indicator for the dependent variable
#'
#' @description
#' Estimates a partial mediation model with multiple categorical indicator for the dependent variable
#'
#' @usage MeasurementYCat(Data,Prior,R)
#'
#' @param Data list(X, M, y_star)
#' @param Prior list(A_M,A_Y)
#' @param R number of MCMC iterations, default = 10000
#'
#' @details
#' *Model*
#'
#' \tabular{ll}{
#' M = beta_0M + Xbeta_1 + U_M  \tab [eq.1] \cr
#' Y = beta_0Y + Mbeta_2 + Xbeta_3 + U_Y \tab [eq.2] \cr
#' }
#'
#' Indicator equations:
#'
#' \tabular{lcl}{
#' y*_1   \tab = \tab M + U_y*_1 \cr
#'  ˜y_1  \tab = \tab  OrdProbit(y*_1,C_y_1) \cr
#'  y*_2    \tab = \tab tau_01 + M + U_y*_2 \cr
#'  ˜y_2  \tab = \tab OrdProbit(y*_2,C_y_2) \cr
#'  ... \tab  \tab  \cr
#'   y*_l   \tab =  \tab tau_0l-1 + M + U_y*_l \cr
#'  ˜y_l  \tab = \tab OrdProbit(y*_l,C_y_l) \cr
#' }
#'
#'
#' *Argument Details*
#'
#' \code{Data = list(X, M, y_star)}
#'
#' \tabular{ll}{
#' \code{X(N x 1) } \tab treatment variable vector \cr
#' \code{M(N x 1) } \tab mediator vector  \cr
#' \code{y_star(N x Y_ind) } \tab dependent variable indicators' matrix \cr
#' }
#'
#' \code{Prior = list(A_M,A_Y)} *[optional]*
#'
#' \tabular{ll}{
#' \code{A_M }   \tab vector of coefficients' prior variances of eq.1, default = rep(100,2) \cr
#' \code{A_Y }   \tab vector of coefficients' prior variances of eq.2, default = c(100,100,1) \cr
#' }
#'
#' @return
#' \tabular{ll}{
#' \code{beta_1(R X 2) } \tab  matrix of eq.1 coefficients' draws \cr
#' \code{beta_2(R X 3) } \tab  matrix of eq.2 coefficients' draws \cr
#' \code{tau (Y_ind X 2 X R) } \tab array of indicator coefficients' draws. Each slice is one draw, where rows represent the indicator equation and columns are the coefficients. All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively. \cr
#' \code{ssq_y_star(R X Y_ind) } \tab  Matrix of indicator equations' coefficients' error variance draws \cr
#' \code{ssq_M(R X 1) } \tab  vector of eq.1 error variance draws \cr
#' \code{mu_draw } \tab  vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor) \cr
#' \code{var_draw } \tab  vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor) \cr
#' }
#' @export
#'
### Description MeasurementYCat estimates a partial mediation model with multiple categorical indicator for the dependent variable
# and observed mediator using a mixture of Metropolis-Hastings and Gibbs sampling
### Arguments:
# Data  list(X, M, y_star)
# Prior list(A_M,A_Y)
# R
### Details:
## Model:
# M = beta_0M + Xbeta_1 + U_M   (eq.1)
# Y = beta_0Y + Mbeta_2 + Xbeta_3 + U_Y  (eq.2)
# indicator equations:
# y*_1 = M + U_y*_1
# ˜y_1 = OrdProbit(y*_1,C_y_1)
# y*_2 = tau_01 + M + U_y*_2
# ˜y_2 = OrdProbit(y*_2,C_y_2)
# ...
# y*_l = tau_0l-1 + M + U_y*_l
# ˜y_l = OrdProbit(y*_l,C_y_l)
## Data = list(X, M, y_star)
# X(N x 1) treatment variable vector
# M(N x 1) mediator vector
# y_star(N x Y_ind) dependent variable indicators' matrix
# Prior = list(A_M,A_Y) [optional]
# A_M vector of coefficients' prior variances of eq.1 (def: rep(100,2))
# A_Y vector of coefficients' prior variances of eq.2 (def: c(100,100,1))
# R number of MCMC iterations (def:10000)
### Value:
# beta_1(R X 2)  matrix of eq.1 coefficients' draws
# beta_2(R X 3)  matrix of eq.2 coefficients' draws
# tau (Y_ind X 2 X R) array of indicator coefficients' draws.
# Each slice is one draw, where rows represent the indicator equation and columns are the coefficients
# All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively
# ssq_y_star(R X Y_ind) Matrix of indicator equations' coefficients' error variance draws
# ssq_M(R X 1) vector of eq.1 error variance draws
# mu_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
# var_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
MeasurementYCat=function(Data,Prior,R=10000){  #,Mcmc){
  #
  # revision history:
  #   3/07  Hsiu-Wen Liu
  #   3/07  fixed naming of dstardraw rossi
  #
  # purpose:
  #   draw from posterior for ordered probit using Gibbs Sampler
  #   and metropolis RW
  #
  # Arguments:
  #   Data - list of X,y,k
  #     X is nobs x nvar, y is nobs vector of 1,2,.,k (ordinal variable)
  #   Prior - list of A, betabar
  #     A is nvar x nvar prior preci matrix
  #     betabar is nvar x 1 prior mean
  #     Ad is ndstar x ndstar prior preci matrix of dstar (ncut is number of cut-offs being estimated)
  #     dstarbar is ndstar x 1 prior mean of dstar
  #   Mcmc
  #     R is number of draws
  #     keep is thinning parameter
  #     nprint - print estimated time remaining on every nprint'th draw
  #     s is scale parameter of random work Metropolis
  #
  # Output:
  #   list of betadraws and cutdraws
  #
  # Model:
  #    z=Xbeta + e  < 0  e ~N(0,1)
  #    y=1,..,k, if z~c(c[k], c[k+1])
  #
  #    cutoffs = c[1],..,c[k+1]
  #    dstar = dstar[1],dstar[k-2]
  #    set c[1]=-100, c[2]=0, ...,c[k+1]=100
  #
  #    c[3]=exp(dstar[1]),c[4]=c[3]+exp(dstar[2]),...,
  #    c[k]=c[k-1]+exp(datsr[k-2])
  #
  # Note: 1. length of dstar = length of cutoffs - 3
  #       2. Be careful in assessing prior parameter, Ad.  .1 is too small for many applications.
  #
  # Prior: beta ~ N(betabar,A^-1)
  #        dstar ~ N(dstarbar, Ad^-1)
  #
  #
  # ----------------------------------------------------------------------
  # Rcpp::sourceCpp('rordprobitGibbs_me_multi_merr.cpp')
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
  if(is.null(Data$y)) {pandterm("Requires Data element y")}
  y=Data$y
  if(is.null(Data$k)) {pandterm("Requires Data element k")}
  k=Data$k
  if(is.null(Data$Y_ind)) {pandterm("Requires Data element Y_ind")}
  Y_ind=Data$Y_ind

  nvar=ncol(X)
  nobs=dim(y)[1]
  ndstar = k-2         # number of dstar being estimated
  ncuts = k+1          # number of cut-offs (including zero and two ends)
  ncut = ncuts-3       # number of cut-offs being estimated c[1]=-100, c[2]=0, c[k+1]=100


  # cutoff_Y_init = Data$cutoff_Y_init
  # Y_tilde_init = Data$Y_tilde_init
  # ssq_y_tilde_init = Data$ssq_y_tilde_init
  # beta_2_init = Data$beta_2_init
  # Y_init = Data$Y_init
  # beta_tilde_init = Data$beta_tilde_init

  #
  # check data for validity
  #
  if(dim(y)[1] != nrow(X) ) {pandterm("y and X not of same row dim")}
  if(  sum(unique(y[,1]) %in% (1:k) ) < length(unique(y[,1])) )   #Tests only y_star1 but I really don't care! :D
  {pandterm("some value of y is not vaild")}

  #
  # check for Prior
  #
  if(missing(Prior))
  { betabar=c(rep(0,nvar)); A=diag(1/c(100,100,1)); A_M=diag(1/rep(100,2)); Ad=diag(ndstar); dstarbar=c(rep(0,ndstar))} #NOTE: A_Y is called A here
  else
  {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))}
    else {betabar=Prior$betabar}
    if(is.null(Prior$A_M)) {A_M=diag(1/rep(100,2))}
    else {A_M=diag(1/Prior$A_M)}
    if(is.null(Prior$A_Y)) {A=diag(1/c(100,100,1))}
    else {A=diag(1/Prior$A_Y)}
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
  #   if(is.null(Mcmc$s)) {s=2.38/sqrt(ndstar)} else {s=Mcmc$s}  #2.38 is the RRscaling
  # }
  #
  # print out problem
  #
  cat(" ", fill=TRUE)
  cat("Starting Gibbs Sampler for Ordered Probit Model",fill=TRUE)
  cat("   with ",nobs,"observations",fill=TRUE)
  cat(" ", fill=TRUE)
  cat("Table of y values",fill=TRUE)
  for(i in 1:Y_ind) print(table(y[,i]))
  cat(" ",fill=TRUE)
  cat("Prior Parms: ",fill=TRUE)
  cat("betabar",fill=TRUE)
  print(betabar)
  cat(" ", fill=TRUE)
  cat("A",fill=TRUE)
  print(A)
  cat(" ", fill=TRUE)
  cat("dstarbar",fill=TRUE)
  print(dstarbar)
  cat(" ", fill=TRUE)
  cat("Ad",fill=TRUE)
  print(Ad)
  cat(" ", fill=TRUE)
  cat("MCMC parms: ",fill=TRUE)
  cat("R= ",R," keep= ",keep," nprint= ",nprint,"s= ",s, fill=TRUE)
  cat(" ",fill=TRUE)

  # use (-Hessian+Ad)^(-1) evaluated at betahat as the basis of the
  # covariance matrix for the random walk Metropolis increments

  betahat = chol2inv(chol(crossprod(X,X)))%*% crossprod(X,rowMeans(y))    ###betahat is computed based on the average of y_stars
  dstarini = c(cumsum(c( rep(0.1, ndstar))))     # set initial value for dstar
  dstarout = optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=X%*%betahat, y=rowMeans(y))
  inc.root=chol(chol2inv(chol((-dstarout$hessian+Ad))))  # chol((H+Ad)^-1)

  ###################################################################
  # Keunwoo Kim
  # 08/20/2014
  ###################################################################
  draws=MeasurementYCatCpp(y,X,k,A,betabar,Ad,s,inc.root,dstarbar,betahat,
                                          Y_ind,
                                          R,keep,nprint)
                                          # cutoff_Y_init, Y_tilde_init, beta_tilde_init, ssq_y_tilde_init, beta_2_init, Y_init)
  ###################################################################
  # draw beta_1, ssq_M | M,X    #This step is independent of the multiple indicator model as we have both M and X
  #
  out<-runiregGibbs_me(Data = list(y=as.matrix(X[,2]),X=cbind(rep(1,nobs),X[,3])),Prior=list(ssq=1,A = A_M),Mcmc = list(R=R))
  beta_1_draw = out$betadraw; ssq_M_draw = out$sigmasqdraw;
  ###################################################################

  draws$cutoff_M=draws$cutoff_M[,2:k,]
  draws$beta_1 = beta_1_draw
  draws$ssq_M = ssq_M_draw
  # attributes(draws$cutdraw)$class="bayesm.mat"
  # attributes(draws$betadraw)$class="bayesm.mat"
  # attributes(draws$dstardraw)$class="bayesm.mat"
  # attributes(draws$cutdraw)$mcpar=c(1,R,keep)
  # attributes(draws$betadraw)$mcpar=c(1,R,keep)
  # attributes(draws$dstardraw)$mcpar=c(1,R,keep)

  return(draws)
}
