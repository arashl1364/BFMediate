#' Sampler for Partial Mediation Model with Multiple Categorical Indicator for the DV
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
#' ## Model
#' \eqn{M = \beta_{0M} + X\beta_1 + U_M}{M = \beta_0M + X\beta_1 + U_M} (eq.1) \cr
#' \eqn{Y = \beta_{0Y} + M\beta_2 + X\beta_3 + U_Y}{Y = \beta_0Y + M\beta_2 + X\beta_3 + U_Y} (eq.2) \cr
#'
#' ## Indicator equations
#' \eqn{y^*_1 = Y + U_{y^*_1}}{y*_1 = Y + U_{y*_1}} \cr
#' \eqn{\tilde{y}_1 = OrdProbit(y^*_1 ,C_{y_1})}{˜y_1 = OrdProbit(y*_1,C_{y_1})} \cr
#' \eqn{y^*_2 = \tau_{20} + Y + U_{y^*_2}}{y*_2 = \tau_20 + Y + U_{y*_2}} \cr
#' \eqn{\tilde{y}_2 = OrdProbit(y^*_2 ,C_{y_2})}{˜y_2 = OrdProbit(y*_2,C_{y_2})} \cr
#' \eqn{...} \cr
#' \eqn{y^*_L = \tau_{L0} + Y + U_{y^*_L}}{y*_L = \tau_L0 + Y + U_{y*_L}} \cr
#' \eqn{\tilde{y}_L = OrdProbit(y^*_L ,C_{y_L})}{˜y_L = OrdProbit(y*_L,C_{y_L})} \cr
#'
#' ## Prior specification:
#'
#'\eqn{\beta_{0M}} \eqn{\sim}{~} N(0,100), \eqn{\beta_{0Y}} \eqn{\sim}{~} \eqn{N(0,100)}  \cr
#'\eqn{\beta_1} \eqn{\sim}{~} \eqn{N(0,100), \beta_2} \eqn{\sim}{~} \eqn{N(0,100), \beta_3} \eqn{\sim}{~} \eqn{N(0,1)}  \cr
#'\eqn{\sigma^2_M} \eqn{\sim}{~} \eqn{Inv\chi^2(\nu,\nu S)} , where \eqn{\nu=1} and S=3. \cr
#'\eqn{\tau_{20}, ..., \tau_{L0}} \eqn{\sim}{~}  \eqn{N(0,100)} \cr
#'\eqn{\sigma^2_{y*_1}, ..., \sigma^2_{y*_L}} \eqn{\sim}{~} \eqn{Inv\chi^2(\nu,\nu S)} \cr
#'\eqn{C*_{y_1}, ..., C*_{y_L}} \eqn{\sim}{~} \eqn{N(0,I)} \cr
#'
#' Note: \eqn{C*_{y_1}, ..., C*_{y_L}} are untransformed
#' cutoffs, which are then exponentially transformed to impose sign and order constraint on them. Subjective prior values for the error variances are \eqn{\nu=1}, S=3.
#'
#' ## Argument Details
#' ## \code{Data = list(X, M, y_star)}
#' \describe{
#'   \item{X(N x 1)}{treatment variable vector}
#'   \item{M(N x 1)}{mediator vector}
#'   \item{y_star(N x Y_ind)}{dependent variable indicators' matrix}
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
#' \item{beta_M(R X 2)}{matrix of eq.1 coefficients' draws }
#' \item{beta_Y(R X 3)}{matrix of eq.2 coefficients' draws }
#' \item{tau(Y_ind X 2 X R)}{array of indicator coefficients' draws. Each slice is one draw, where rows represent the indicator equation and columns the coefficients. All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively. }
#' \item{ssq_y_star(R X Y_ind)}{  Matrix of indicator equations' coefficients' error variance draws }
#' \item{ssq_M(R X 1)}{vector of eq.1 error variance draws }
#' \item{cutoff_Y (Y_ind X k_Y X R)}{array of discretized dependent variable indicators' cutoff values.}
#' \item{Ydraw(R X N)}{Matrix of the augmented latent dependent variable}
#' \item{mu_draw}{vector of means indexing MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#' \item{var_draw}{vector of variance indexing MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#' }
#' @export
#' @examples
#' set.seed(60)
#' SimMeasurementYCat = function(X, beta_M, beta_Y, sigma_M, cutoff_Y, Y_ind, tau, ssq_y_star){
#'
#'   nobs = dim(X)[1]
#'   y_tilde = y_star = matrix(double(nobs*Y_ind), ncol = Y_ind)
#'
#'   M = beta_M[1] + beta_M[2] * X + rnorm(nobs) * sigma_M
#'   Y = beta_Y[1] + beta_Y[2] * M + beta_Y[3] * X + rnorm(nobs)
#'
#'   for(i in 1: Y_ind){
#'     y_star[,i] = tau[i] + Y + sqrt(ssq_y_star[i])*rnorm(nobs);
#'     y_tilde[,i] = cut(y_star[,i], br = cutoff_Y[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
#'   }
#'
#'   return(list(Y = Y, M = M, y_tilde = y_tilde, X = X,
#'               beta_M = beta_M,
#'               k_Y=dim(cutoff_Y)[2]-1, beta_Y = beta_Y, tau = tau,
#'              ssq_y_star = ssq_y_star, y_star = y_star, cutoff_Y = cutoff_Y,
#'               Y_ind=Y_ind))
#' }
#'
#' Y_ind = 2
#' Ycut = 8
#' nobs = 1000
#' X=as.matrix(runif(nobs,min=0, max=1))
#' beta_M = c(.5,1)
#' beta_Y = c(1, 4, 2)
#' sigma_M = 1^.5
#' ssq_y_star = c(.5,.7)
#' tau = c(0,-.5)   #first intercept should always be 0
#' cutoff_Y =  matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
#'                      -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Ycut, byrow = TRUE)
#' DataYCat = SimMeasurementYCat(X, beta_M, beta_Y, sigma_M, cutoff_Y, Y_ind, tau, ssq_y_star)
#'
#' #estimation
#' A_M = c(100,100); #Prior variance for beta_0M, beta_1
#' A_Y = c(100,100,1) #Prior variance for beta_0Y, beta_2, beta_3
#' Prior = list(A_M = A_M, A_Y = A_Y)
#' Ycut = max(as.matrix(DataYCat$y_tilde)[,1]) +1
#' Data = list(X=cbind(rep(1,length(DataYCat$X)),DataYCat$M,DataYCat$X), y = as.matrix(DataYCat$y_tilde),
#'             k=Ycut-1, Y_ind=dim(as.matrix(DataYCat$y_tilde))[2])
#' out = MeasurementYCat(Data=Data, Prior=Prior, R=10000)
#'
#' #results
#' colMeans(out$beta_M)
#' colMeans(out$beta_Y)
#' apply(out$cutoff_Y,c(1,2),FUN = mean)
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
# beta_M(R X 2)  matrix of eq.1 coefficients' draws
# beta_Y(R X 3)  matrix of eq.2 coefficients' draws
# tau (Y_ind X 2 X R) array of indicator coefficients' draws.
# Each slice is one draw, where rows represent the indicator equation and columns are the coefficients
# All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively
# ssq_y_star(R X Y_ind) Matrix of indicator equations' coefficients' error variance draws
# ssq_M(R X 1) vector of eq.1 error variance draws
# mu_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
# var_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
MeasurementYCat=function(Data,Prior,R=10000){  #,Mcmc){
  #
  # A. Laghaie 2019
  # The ordered probit part of the inference and the R shell is based on rordprobitGibbs {bayesm}
  #
  # ----------------------------------------------------------------------
  # define functions needed
  #  dstartoc is a fuction to transfer dstar to its cut-off value

  dstartoc=function(dstar) {c(-100, 0, cumsum(exp(dstar)), 100)}

  # compute conditional likelihood of data given cut-offs
  #
  lldstar=function(dstar,y,mu){
    gamma=dstartoc(dstar)
    arg = stats::pnorm(gamma[y+1]-mu)-stats::pnorm(gamma[y]-mu)
    epsilon=1.0e-50
    arg=ifelse(arg < epsilon,epsilon,arg)
    return(sum(log(arg)))
  }
  #
  # ----------------------------------------------------------------------
  #
  # check arguments
  #
  if(missing(Data)) {stop("Requires Data argument -- list of y and X")}
  if(is.null(Data$X)) {stop("Requires Data element X")}
  X=Data$X
  if(is.null(Data$y)) {stop("Requires Data element y")}
  y=Data$y
  if(is.null(Data$k)) {stop("Requires Data element k")}
  k=Data$k
  if(is.null(Data$Y_ind)) {stop("Requires Data element Y_ind")}
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
  if(dim(y)[1] != nrow(X) ) {stop("y and X not of same row dim")}
  if(  sum(unique(y[,1]) %in% (1:k) ) < length(unique(y[,1])) )   #Tests only y_star1 but I really don't care! :D
  {stop("some value of y is not vaild")}

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
  {stop(paste("bad dimensions for A",dim(A)))}
  if(length(betabar) != nvar)
  {stop(paste("betabar wrong length, length= ",length(betabar)))}
  if(ncol(Ad) != nrow(Ad) || ncol(Ad) != ndstar || nrow(Ad) != ndstar)
  {stop(paste("bad dimensions for Ad",dim(Ad)))}
  if(length(dstarbar) != ndstar)
  {stop(paste("dstarbar wrong length, length= ",length(dstarbar)))}

  #
  # check MCMC argument
  #
  keep = 1
  nprint = 100
  s = 2.38/sqrt(ndstar)
  # if(missing(Mcmc)) {stop("requires Mcmc argument")}
  # else
  # {
  #   if(is.null(Mcmc$R))
  #   {stop("requires Mcmc element R")} else {R=Mcmc$R}
  #   if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
  #   if(is.null(Mcmc$nprint)) {nprint=100} else {nprint=Mcmc$nprint}
  #   if(nprint<0) {stop('nprint must be an integer greater than or equal to 0')}
  #   if(is.null(Mcmc$s)) {s=2.38/sqrt(ndstar)} else {s=Mcmc$s}  #2.38 is the RRscaling
  # }
  #
  # print out problem
  #
  cat(" ", fill=TRUE)
  cat("Starting Gibbs Sampler for mediation LVM Model",fill=TRUE)
  cat("   with ",nobs,"observations",fill=TRUE)
  cat(" ", fill=TRUE)
  cat("Table of y values",fill=TRUE)
  for(i in 1:Y_ind) print(table(y[,i]))
  cat(" ",fill=TRUE)
  # cat("Prior Parms: ",fill=TRUE)
  # cat("betabar",fill=TRUE)
  # print(betabar)
  # cat(" ", fill=TRUE)
  # cat("A",fill=TRUE)
  # print(A)
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
  dstarout = stats::optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=X%*%betahat, y=rowMeans(y))
  inc.root=chol(chol2inv(chol((-dstarout$hessian+Ad))))  # chol((H+Ad)^-1)

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

  draws$cutoff_Y=draws$cutoff_Y[,2:k,]
  draws$beta_M = beta_1_draw
  draws$ssq_M = ssq_M_draw

  return(draws)
}
