#' asdfsaf
#'
#' @param Data asdf
#' @param Prior asdf
#' @param Mcmc aasdf
#'
#' @return
#' @export
#'
#' @import bayesm
Mediation_Ordered_Multi_Merr=function(Data,Prior,Mcmc){
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
  #   Data - list of X,y_star,k_Y
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
  lldstar=function(dstar_Y,y_star,mu){
    gamma=dstartoc(dstar_Y)
    arg = pnorm(gamma[y_star+1]-mu)-pnorm(gamma[y_star]-mu)
    epsilon=1.0e-50
    arg=ifelse(arg < epsilon,epsilon,arg)
    return(sum(log(arg)))
  }
  #
  # ----------------------------------------------------------------------
  #
  # check arguments
  #
  if(missing(Data)) {pandterm("Requires Data argument -- list of y_star and X")}
  if(is.null(Data$X)) {pandterm("Requires Data element X")}
  X=Data$X
  if(is.null(Data$y_star)) {pandterm("Requires Data element y_star")}
  y_star=Data$y_star
  if(is.null(Data$k_Y)) {pandterm("Requires Data element k_Y")}
  k_Y=Data$k_Y
  if(is.null(Data$m_star)) {pandterm("Requires Data element m_star")}
  m_star=Data$m_star
  if(is.null(Data$k_M)) {pandterm("Requires Data element k_M")}
  k_M=Data$k_M
  if(is.null(Data$M_ind)) {pandterm("Requires Data element M_ind")}
  M_ind=Data$M_ind
  if(is.null(Data$Y_ind)) {pandterm("Requires Data element Y_ind")}
  Y_ind=Data$Y_ind

  #fixing parameters for testing the sampler
  cutoffs_M_init = Data$cutoffs_M_init
  M_tilde_init = Data$M_tilde_init
  beta_m_tilde_init = Data$beta_m_tilde_init
  ssq_m_tilde_init = Data$ssq_m_tilde_init
  beta_init = Data$beta_init
  M_init = Data$M_init

  cutoffs_Y_init = Data$cutoffs_Y_init
  Y_tilde_init = Data$Y_tilde_init
  beta_y_tilde_init = Data$beta_y_tilde_init
  ssq_y_tilde_init = Data$ssq_y_tilde_init
  beta_2_init = Data$beta_2_init
  Y_init = Data$Y_init

  # cutoffs_M_init,
  # M_tilde_init, beta_tilde_init, ssq_m_tilde_init, beta_init, M_init,
  # cutoffs_Y_init,
  # Y_tilde_init, beta_2_tilde_init, ssq_y_tilde_init, beta_2_init, Y_init)



  nvar_Y=ncol(X)+1
  nobs=dim(y_star)[1]
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
  if(dim(y_star)[1] != nrow(X) ) {pandterm("y_star and X not of same row dim")}
  if(  sum(unique(y_star[,1]) %in% (1:k_Y) ) < length(unique(y_star[,1])) )     #Tests only y_star1 but I really don't care! :D
  {pandterm("some value of y_star is not vaild")}
  if(dim(m_star)[1] != nrow(X) ) {pandterm("m_star and X not of same row dim")}
  if(  sum(unique(m_star[,1]) %in% (1:k_M) ) < length(unique(m_star[,1])) )     #Tests only m_star1 but I really don't care! :D
  {pandterm("some value of m_star is not vaild")}

  #
  # check for Prior
  #
  if(missing(Prior))
  { betabar=c(rep(0,nvar_M)); A_M=.01*diag(nvar_M); Ad_M=diag(ndstar_M); dstarbar_M=c(rep(0,ndstar_M))
  beta_2_bar=c(rep(0,nvar_Y)); A_Y=.01*diag(nvar_Y); Ad_Y=diag(ndstar_Y); dstarbar_Y=c(rep(0,ndstar_Y))}
  else
  {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar_M))}
    else {betabar=Prior$betabar}
    if(is.null(Prior$A_M)) {A_M=.01*diag(nvar_M)}
    else {A_M=Prior$A_M}
    if(is.null(Prior$Ad_M)) {Ad_M=diag(ndstar_M)}
    else {Ad_M=Prior$Ad_M}
    if(is.null(Prior$dstarbar_M)) {dstarbar_M=c(rep(0,ndstar_M))}
    else {dstarbar_M=Prior$dstarbar_M}

    if(is.null(Prior$beta_2_bar)) {beta_2_bar=c(rep(0,nvar_Y))}
    else {beta_2_bar=Prior$beta_2_bar}
    if(is.null(Prior$A_Y)) {A_Y=.01*diag(nvar_Y)}
    else {A_Y=Prior$A_Y}
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
  if(missing(Mcmc)) {pandterm("requires Mcmc argument")}
  else
  {
    if(is.null(Mcmc$R))
    {pandterm("requires Mcmc element R")} else {R=Mcmc$R}
    if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
    if(is.null(Mcmc$nprint)) {nprint=100} else {nprint=Mcmc$nprint}
    if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
    if(is.null(Mcmc$s_M)) {s_M=2.38/sqrt(ndstar_M)} else {s_M=Mcmc$s_M} #2.38 is RRscaling
    if(is.null(Mcmc$s_Y)) {s_Y=2.38/sqrt(ndstar_Y)} else {s_Y=Mcmc$s_Y}
  }
  #
  # print out problem
  #
  cat(" ", fill=TRUE)
  cat("Starting Gibbs Sampler for Ordered Probit Model",fill=TRUE)
  cat("   with ",nobs,"observations",fill=TRUE)
  cat(" ", fill=TRUE)
  cat("Table of y_star values",fill=TRUE)
  for(i in 1:Y_ind) print(table(y_star[,i]))
  cat("Table of m_star values",fill=TRUE)
  for(i in 1:M_ind) print(table(m_star[,i]))
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

  betahat = chol2inv(chol(crossprod(X,X)))%*% crossprod(X,rowMeans(m_star))
  dstarini = c(cumsum(c( rep(0.1, ndstar_M))))     # set initial value for dstar_M
  dstarout = optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=X%*%betahat, y=rowMeans(m_star))
  inc.root_M=chol(chol2inv(chol((-dstarout$hessian+Ad_M))))  # chol((H+Ad_Y)^-1)


  # use (-Hessian+Ad_Y)^(-1) evaluated at beta_2_hat as the basis of the
  # covariance matrix for the random walk Metropolis increments

  Xmstar=cbind(X[,1],rowMeans(m_star),X[,2])
  beta_2_hat = chol2inv(chol(crossprod(Xmstar,Xmstar)))%*% crossprod(Xmstar,rowMeans(y_star))
  dstarini = c(cumsum(c( rep(0.1, ndstar_Y))))     # set initial value for dstar_Y
  dstarout = optim(dstarini, lldstar, method = "BFGS", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=Xmstar%*%beta_2_hat, y=rowMeans(y_star))
  inc.root_Y=chol(chol2inv(chol((-dstarout$hessian+Ad_Y))))  # chol((H+Ad_Y)^-1)

  ###################################################################
  ###################################################################
  # Keunwoo Kim
  # 08/20/2014
  # Modified by Arash Laghaie
  # 12/27/2018
  ###################################################################
  draws= Mediation_Ordered_Multi_Merr_cpp(X, m_star, y_star, k_M, k_Y,M_ind,Y_ind,
                                     A_M, betabar, Ad_M, s_M, inc.root_M, dstarbar_M, betahat,
                                     A_Y, beta_2_bar, Ad_Y, s_Y, inc.root_Y, dstarbar_Y, beta_2_hat,
                                     R, keep, nprint)
                                     # cutoffs_M_init,
                                     # M_tilde_init, beta_m_tilde_init, ssq_m_tilde_init, beta_init, M_init,
                                     # cutoffs_Y_init,
                                     # Y_tilde_init, beta_y_tilde_init, ssq_y_tilde_init, beta_2_init, Y_init)

  ###################################################################

  draws$cutdraw_M=draws$cutdraw_M[,2:k_M,]
  # attributes(draws$cutdraw_M)$class="bayesm.mat"
  # attributes(draws$betadraw)$class="bayesm.mat"
  # attributes(draws$dstardraw_M)$class="bayesm.mat"
  # attributes(draws$cutdraw_M)$mcpar=c(1,R,keep)
  # attributes(draws$betadraw)$mcpar=c(1,R,keep)
  # attributes(draws$dstardraw_M)$mcpar=c(1,R,keep)

  draws$cutdraw_Y=draws$cutdraw_Y[,2:k_Y,]
  # attributes(draws$cutdraw_Y)$class="bayesm.mat"
  # attributes(draws$beta_2_draw)$class="bayesm.mat"
  # attributes(draws$dstardraw_Y)$class="bayesm.mat"
  # attributes(draws$cutdraw_Y)$mcpar=c(1,R,keep)
  # attributes(draws$beta_2_draw)$mcpar=c(1,R,keep)
  # attributes(draws$dstardraw_Y)$mcpar=c(1,R,keep)

  return(draws)
}
