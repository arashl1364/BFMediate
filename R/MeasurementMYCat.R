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
#' ## Model
#' \eqn{M = \beta_{0M} + X\beta_1 + U_M}{M = \beta_0M + X\beta_1 + U_M} (eq.1) \cr
#' \eqn{Y = \beta_{0Y} + M\beta_2 + X\beta_3 + U_Y}{Y = \beta_0Y + M\beta_2 + X\beta_3 + U_Y} (eq.2) \cr
#'
#' ## Measurement equations
#' \eqn{m^*_1 = M + U_{m^*_1}}{m*_1 = M + U_{m*_1}} \cr
#' \eqn{\tilde{m}_1 = OrdProbit(m^*_1 ,C_{m_1})}{˜m_1 = OrdProbit(m*_1,C_{m_1})} \cr
#' \eqn{m^*_2 = \lambda_{20} + M + U_{m^*_2}}{m*_2 = \lambda_20 + M + U_{m*_2}} \cr
#' \eqn{\tilde{m}_2 = OrdProbit(m^*_2 ,C_{m_2})}{˜m_2 = OrdProbit(m*_2,C_{m_2})} \cr
#' \eqn{...} \cr
#' \eqn{m^*_K = \lambda_{K0} + M + U_{m^*_K}}{m*_K = \lambda_K0 + M + U_{m*_K}} \cr
#' \eqn{\tilde{m}_K = OrdProbit(m^*_K ,C_{m_K})}{˜m_K = OrdProbit(m*_K,C_{m_K})} \cr
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
#'\eqn{\lambda_{20},...,\lambda_{K0}} \eqn{\sim}{~}  \eqn{N(0,100)}  \cr
#'\eqn{\sigma^2_{m*_1}, ..., \sigma^2_{m*_K}} \eqn{\sim}{~} \eqn{Inv\chi^2(\nu,\nu S)} \cr
#'\eqn{\tau_{20}, ..., \tau_{L0}} \eqn{\sim}{~}  \eqn{N(0,100)} \cr
#'\eqn{\sigma^2_{y*_1}, ..., \sigma^2_{y*_L}} \eqn{\sim}{~} \eqn{Inv\chi^2(\nu,\nu S)} \cr
#'\eqn{C*_{m_1}, ..., C*_{m_K}, C*_{y_1}, ..., C*_{y_L}} \eqn{\sim}{~} \eqn{N(0,I)} \cr
#'
#' Note: \eqn{C*_{m_1}, ..., C*_{m_K}, C*_{y_1}, ..., C*_{y_L}} are untransformed
#' cutoffs, which are then exponentially transformed to impose sign and order constraint on them. Subjective prior values for the error variances are \eqn{\nu=1}, S=3.
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
#'   \item{beta_M(R X 2)}{matrix of eq.1 coefficients' draws}
#'   \item{beta_Y(R X 3)}{matrix of eq.2 coefficients' draws}
#'   \item{lambda (M_ind X 2 X R)}{array of mediator indicator coefficients' draws. Each slice is one draw, where rows represent the indicator equation and columns the coefficients. All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively.}
#'   \item{tau (Y_ind X 2 X R)}{array of dependent variable indicator coefficients' draws. Each slice is one draw, where rows represent the indicator equation and columns the coefficients. All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively.}
#'   \item{ssq_m_star(R X M_ind)}{Matrix of mediator indicator equations' coefficients' error variance draws}
#'   \item{ssq_y_star(R X Y_ind)}{Matrix of dependent variable indicator equations' coefficients' error variance draws}
#'   \item{cutoff_M (M_ind X k_M X R)}{array of discretized mediator indicators' cutoff values.}
#'   \item{cutoff_Y (Y_ind X k_Y X R)}{array of discretized dependent variable indicators' cutoff values.}
#'   \item{Mdraw(R X N)}{Matrix of the augmented latent mediator}
#'   \item{Ydraw(R X N)}{Matrix of the augmented latent dependent variable}
#'   \item{mu_draw}{vector of means indexing MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#'   \item{var_draw}{vector of variance indexing MCMC draws of the direct effect (used in BFSD to compute Bayes factor)}
#' }
#'
#' @export
#' @examples
#' set.seed(60)
#' SimMeasurementMYCat = function(X, beta_M, cutoff_M, beta_Y, cutoff_Y, M_ind, Y_ind,
#'                                lambda, tau, ssq_m_star, ssq_y_star){
#'
#'   nobs = dim(X)[1]
#'   m_tilde = m_star = matrix(double(nobs*M_ind), ncol = M_ind)
#'   y_tilde = y_star = matrix(double(nobs*Y_ind), ncol = Y_ind)
#'
#'   M = cbind(rep(1,nobs),X)%*%beta_M + rnorm(nobs)
#'   for(i in 1: M_ind){
#'     m_star[,i] = lambda[i] + M + sqrt(ssq_m_star[i])*rnorm(nobs);
#'     m_tilde[,i] = cut(m_star[,i], br = cutoff_M[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
#'   }
#'
#'   Y = cbind(rep(1,nobs),cbind(M,X))%*%beta_Y + rnorm(nobs)
#'   for(i in 1: Y_ind){
#'     y_star[,i] = tau[i] + Y + sqrt(ssq_y_star[i])*rnorm(nobs);
#'     y_tilde[,i] = cut(y_star[,i], br = cutoff_Y[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
#'   }
#'
#'   return(list(Y = Y, M = M, y_tilde = y_tilde, m_tilde = m_tilde, X = X,
#'               k_M=dim(cutoff_M)[2]-1, beta_M = beta_M, lambda = lambda,
#'               ssq_m_star = ssq_m_star, m_star = m_star, cutoff_M = cutoff_M,
#'               k_Y=dim(cutoff_Y)[2]-1, beta_Y = beta_Y, tau = tau, ssq_y_star = ssq_y_star,
#'               y_star = y_star, cutoff_Y = cutoff_Y, M_ind=M_ind, Y_ind=Y_ind))
#' }
#'
#' M_ind = 2
#' Y_ind = 2
#' Mcut = Ycut = 8
#' nobs=1000
#' X=as.matrix(runif(nobs,min=0, max=1))
#' beta_M = c(.5,1)
#' beta_Y = c(.7, 1.5, 0)
#' ssq_m_star = c(.5,.3)
#' lambda = c(0,-.5) #the intercepts for the latent M indicators w. measurement error
#' ssq_y_star = c(.2,.2)
#' tau = c(0,-.5)
#' cutoff_M = matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
#'                     -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Mcut, byrow = TRUE)
#' cutoff_Y =  matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
#'                      -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= Ycut, byrow = TRUE)
#' DataMYCat = SimMeasurementMYCat(X, beta_M, cutoff_M, beta_Y, cutoff_Y, M_ind,
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
#' colMeans(out$beta_M)
#' colMeans(out$beta_Y)
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
# beta_M(R X 2)  matrix of eq.1 coefficients' draws
# beta_Y(R X 3)  matrix of eq.2 coefficients' draws
# lambda (M_ind X 2 X R) array of mediator indicator coefficients' draws.
# tau (Y_ind X 2 X R) array of dependent variable indicator coefficients' draws.
# Each slice is one draw, where rows represent the indicator equation and columns the coefficients
# All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively
# ssq_m_star(R X M_ind) Matrix of mediator indicator equations' coefficients' error variance draws
# ssq_y_star(R X Y_ind) Matrix of dependent variable indicator equations' coefficients' error variance draws
# mu_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
# var_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)

MeasurementMYCat=function(Data,Prior,R=10000){
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
  lldstar=function(dstar_Y,y_tilde,mu){
    gamma=dstartoc(dstar_Y)
    arg = stats::pnorm(gamma[y_tilde+1]-mu)-stats::pnorm(gamma[y_tilde]-mu)
    epsilon=1.0e-50
    arg=ifelse(arg < epsilon,epsilon,arg)
    return(sum(log(arg)))
  }
  #
  # ----------------------------------------------------------------------
  #
  # check arguments
  #
  if(missing(Data)) {stop("Requires Data argument -- list of y_tilde and X")}
  if(is.null(Data$X)) {stop("Requires Data element X")}
  X=Data$X
  if(is.null(Data$y_tilde)) {stop("Requires Data element y_tilde")}
  y_tilde=Data$y_tilde
  if(is.null(Data$k_Y)) {stop("Requires Data element k_Y")}
  k_Y=Data$k_Y
  if(is.null(Data$m_tilde)) {stop("Requires Data element m_tilde")}
  m_tilde=Data$m_tilde
  if(is.null(Data$k_M)) {stop("Requires Data element k_M")}
  k_M=Data$k_M
  if(is.null(Data$M_ind)) {stop("Requires Data element M_ind")}
  M_ind=Data$M_ind
  if(is.null(Data$Y_ind)) {stop("Requires Data element Y_ind")}
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
  if(dim(y_tilde)[1] != nrow(X) ) {stop("y_tilde and X not of same row dim")}
  if(  sum(unique(y_tilde[,1]) %in% (1:k_Y) ) < length(unique(y_tilde[,1])) )     #Tests only y_star1 but I really don't care! :D
  {stop("some value of y_tilde is not vaild")}
  if(dim(m_tilde)[1] != nrow(X) ) {stop("m_tilde and X not of same row dim")}
  if(  sum(unique(m_tilde[,1]) %in% (1:k_M) ) < length(unique(m_tilde[,1])) )     #Tests only m_star1 but I really don't care! :D
  {stop("some value of m_tilde is not vaild")}

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
  {stop(paste("bad dimensions for A_M",dim(A_M)))}
  if(length(betabar) != nvar_M)
  {stop(paste("betabar wrong length, length= ",length(betabar)))}
  if(ncol(Ad_M) != nrow(Ad_M) || ncol(Ad_M) != ndstar_M || nrow(Ad_M) != ndstar_M)
  {stop(paste("bad dimensions for Ad_M",dim(Ad_M)))}
  if(length(dstarbar_M) != ndstar_M)
  {stop(paste("dstarbar_M wrong length, length= ",length(dstarbar_M)))}

  if(ncol(A_Y) != nrow(A_Y) || ncol(A_Y) != nvar_Y || nrow(A_Y) != nvar_Y)
  {stop(paste("bad dimensions for A_Y",dim(A_Y)))}
  if(length(beta_2_bar) != nvar_Y)
  {stop(paste("beta_2_bar wrong length, length= ",length(beta_2_bar)))}
  if(ncol(Ad_Y) != nrow(Ad_Y) || ncol(Ad_Y) != ndstar_Y || nrow(Ad_Y) != ndstar_Y)
  {stop(paste("bad dimensions for Ad_Y",dim(Ad_Y)))}
  if(length(dstarbar_Y) != ndstar_Y)
  {stop(paste("dstarbar_Y wrong length, length= ",length(dstarbar_Y)))}

  #
  # check MCMC argument
  #
  keep = 1
  nprint = 100
  s_M = 2.38/sqrt(ndstar_M)
  s_Y = 2.38/sqrt(ndstar_Y)
  # if(missing(Mcmc)) {stop("requires Mcmc argument")}
  # else
  # {
  #   if(is.null(Mcmc$R))
  #   {stop("requires Mcmc element R")} else {R=Mcmc$R}
  #   if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
  #   if(is.null(Mcmc$nprint)) {nprint=100} else {nprint=Mcmc$nprint}
  #   if(nprint<0) {stop('nprint must be an integer greater than or equal to 0')}
  #   if(is.null(Mcmc$s_M)) {s_M=2.38/sqrt(ndstar_M)} else {s_M=Mcmc$s_M} #2.38 is RRscaling
  #   if(is.null(Mcmc$s_Y)) {s_Y=2.38/sqrt(ndstar_Y)} else {s_Y=Mcmc$s_Y}
  # }
  #
  # print out problem
  #
  cat(" ", fill=TRUE)
  cat("Starting Gibbs Sampler for mediation LVM Model",fill=TRUE)
  cat("   with ",nobs,"observations",fill=TRUE)
  cat(" ", fill=TRUE)
  cat("Table of y_tilde values",fill=TRUE)
  for(i in 1:Y_ind) print(table(y_tilde[,i]))
  cat(" ", fill=TRUE)
  cat("Table of m_tilde values",fill=TRUE)
  for(i in 1:M_ind) print(table(m_tilde[,i]))
  cat(" ",fill=TRUE)
  # cat("Prior Parms: ",fill=TRUE)
  # cat("betabar",fill=TRUE)
  # print(betabar)
  # cat("beta_2_bar",fill=TRUE)
  # print(beta_2_bar)
  # cat(" ", fill=TRUE)
  # cat("A_M",fill=TRUE)
  # print(A_M)
  # cat("A_Y",fill=TRUE)
  # print(A_Y)
  # cat(" ", fill=TRUE)
  # cat("dstarbar_M",fill=TRUE)
  # print(dstarbar_M)
  # cat("dstarbar_Y",fill=TRUE)
  # print(dstarbar_Y)
  # cat(" ", fill=TRUE)
  # cat("Ad_M",fill=TRUE)
  # print(Ad_M)
  # cat("Ad_Y",fill=TRUE)
  # print(Ad_Y)
  # cat(" ", fill=TRUE)
  # cat("MCMC parms: ",fill=TRUE)
  # cat("R= ",R," keep= ",keep," nprint= ",nprint,"s_M= ",s_M,"s_Y= ",s_Y, fill=TRUE)
  # cat(" ",fill=TRUE)

  # use (-Hessian+Ad_M)^(-1) evaluated at betahat as the basis of the
  # covariance matrix for the random walk Metropolis increments

  betahat = chol2inv(chol(crossprod(X,X)))%*% crossprod(X,rowMeans(m_tilde))
  dstarini = c(cumsum(c( rep(0.1, ndstar_M))))     # set initial value for dstar_M
  dstarout = stats::optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=X%*%betahat, y=rowMeans(m_tilde))
  inc.root_M=chol(chol2inv(chol((-dstarout$hessian+Ad_M))))  # chol((H+Ad_Y)^-1)


  # use (-Hessian+Ad_Y)^(-1) evaluated at beta_2_hat as the basis of the
  # covariance matrix for the random walk Metropolis increments

  Xmstar=cbind(X[,1],rowMeans(m_tilde),X[,2])
  beta_2_hat = chol2inv(chol(crossprod(Xmstar,Xmstar)))%*% crossprod(Xmstar,rowMeans(y_tilde))
  dstarini = c(cumsum(c( rep(0.1, ndstar_Y))))     # set initial value for dstar_Y
  dstarout = stats::optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=Xmstar%*%beta_2_hat, y=rowMeans(y_tilde))
  inc.root_Y=chol(chol2inv(chol((-dstarout$hessian+Ad_Y))))  # chol((H+Ad_Y)^-1)


  ###################################################################
  draws= MeasurementMYCatCpp(X, m_tilde, y_tilde, k_M, k_Y,M_ind,Y_ind,
                                     A_M, betabar, Ad_M, s_M, inc.root_M, dstarbar_M, betahat,
                                     A_Y, beta_2_bar, Ad_Y, s_Y, inc.root_Y, dstarbar_Y, beta_2_hat,
                                     R, keep, nprint)
  ###################################################################

  draws$cutoff_M=draws$cutoff_M[,2:k_M,]

  draws$cutoff_Y=draws$cutoff_Y[,2:k_Y,]

  return(draws)
}
