#' Partial Mediation
#'
#' @param Data Data
#' @param pars Parameters
#' @param R Dunno what is it
#'
#' @return something
#' @export
#'
#Description PartialMed estimates a partial mediation model using series of Gibbs Samplers
#Arguments:
# Data  list(X, M, Y)
# pars list(A_M,A_Y)
# R
#Details:
# Data = list(X, M, Y)
# X(N x 1) treatment variable vector
# M(N x 1) mediator vector
# Y(N x 1) dependent variable vector
# pars = list(A_M,A_Y) [optional]
# A_M vector of coefficients' prior variances of eq.1 (def: rep(100,2))
# A_Y vector of coefficients' prior variances of eq.2 (def: c(100,100,1))
# R number of MCMC iterations (def:10000)
#Value:
# beta_1(R X 2)  matrix of eq.1 coefficients' posterior draws
# beta_2(R X 3)  matrix of eq.2 coefficients' posterior draws
# ssq_M(R X 1) vector of eq.1 error variance posterior draws
# ssq_Y(R X 1) vector of eq.2 error variance posterior draws
# mu_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
# var_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)

PartialMed=function(Data, pars, R=10000){

    #initialization and memory allocation

  #Data
  X=as.matrix(Data$X); Y=Data$Y; M=Data$M;
  N=length(Y)
  k=dim(X)[2] + 1  # accounting for the intercept
  if(is.null(k)) k=2

  if(missing(pars)){
    A_M = 1/rep(100,2)  #rep(.01,2)
    A_Y = 1/c(100,100,1) #rep(.01,3)
  }
  else
  {
    if(is.null(pars$A_M)) {A_M = 1/rep(100,2)} #{A_M = rep(.01,2)}
    else {A_M = 1/pars$A_M}
    if(is.null(pars$A_Y)) {A_Y = 1/c(100,100,1)} #{A_Y = rep(.01,3)}
    else {A_Y = 1/pars$A_Y}
  }

  ##Posterior draws
  ssq_M_draw=ssq_y_draw=rep(0,R)
  beta_1_draw=matrix(double(R*k),ncol = k)
  beta_2_draw=matrix(double(R*(k+1)),ncol = k+1) # beta_2 is c(intercept, beta2, beta3)

  ##Moment draws
  mu_beta_1_draw = matrix(double(R*k),ncol = k)
  mu_beta_2_draw = matrix(double(R*(k+1)),ncol = k+1)
  IR_beta_1_draw = array(double((k^2)*R), dim = c(k,k,R))
  IR_beta_2_draw =  array(double(((k+1)^2)*R), dim = c(k+1,k+1,R))
  nu_ssq_M_draw = S_ssq_M_draw = nu_ssq_y_draw = S_ssq_y_draw = rep(0,R)

  #ssq_M=ssq_y=1;

  itime=proc.time()[3]
  cat("MCMC Iteration (est time to end - min) ",fill=TRUE)


    # draw beta_1, ssq_M | M,X
    #
    out<-runiregGibbs_me(Data = list(y=M,X=cbind(rep(1,N),X)),Prior=list(ssq=1,A = as.matrix(diag(A_M,k))),Mcmc = list(R=R))
    beta_1_draw = out$betadraw; ssq_M_draw = out$sigmasqdraw;
    # ##Moments
    # mu_beta_1_draw = out$mubeta
    # IR_beta_1_draw = out$IR      #covariance matrix of beta_1
    # nu_ssq_M_draw = out$nu
    # S_ssq_M_draw = out$S
    #
    #
    # draw beta_2, ssq_y | Y,M,X
    #
    out<-runiregGibbs_me(Data = list(y=Y,X=cbind(rep(1,N),M,X)),Prior=list(ssq=1, A =  as.matrix(diag(A_Y,k+1))), Mcmc = list(R=R))
    beta_2_draw = out$betadraw; ssq_y_draw = out$sigmasqdraw;
    ##Moments
    mu_beta_2_draw = out$mubeta
    IR_beta_2_draw = out$IR    #covariance matrix of beta draws
    # nu_ssq_y_draw = out$nu
    # S_ssq_y_draw = out$S

  ctime = proc.time()[3]
  cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

  # attributes(beta_1_draw)$class=c("bayesm.mat","mcmc")
  # attributes(beta_1_draw)$mcpar=c(1,R,1)
  # attributes(beta_2_draw)$class=c("bayesm.mat","mcmc")
  # attributes(beta_2_draw)$mcpar=c(1,R,1)
  # attributes(ssq_M_draw)$class=c("bayesm.var","bayesm.mat","mcmc")
  # attributes(ssq_M_draw)$mcpar=c(1,R,1)
  # attributes(ssq_y_draw)$class=c("bayesm.var","bayesm.mat","mcmc")
  # attributes(ssq_y_draw)$mcpar=c(1,R,1)

  return(list(beta_1=beta_1_draw,beta_2=beta_2_draw, ssq_M=ssq_M_draw ,ssq_y=ssq_y_draw,
              mu_draw=mu_beta_2_draw[,3], var_draw=IR_beta_2_draw[3,3,]))    #MCMC moments of the direct effect

  # return(list(beta_1=beta_1_draw,beta_2=beta_2_draw, ssq_M=ssq_M_draw ,ssq_y=ssq_y_draw,
  #             mu_beta_1=mu_beta_1_draw, IR_beta_1=IR_beta_1_draw, nu_ssq_M=nu_ssq_M_draw,S_ssq_M=S_ssq_M_draw,
  #             mu_beta_2=mu_beta_2_draw, IR_beta_2=IR_beta_2_draw, nu_ssq_y=nu_ssq_y_draw, S_ssq_y=S_ssq_y_draw))
}


