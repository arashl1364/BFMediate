#' MeasurementCont
#'
#' @param Data asd
#' @param Prior asd
#' @param Mcmc asd
#'
#' @return asd
#' @export
#'
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
# beta_1(R X 2)  matrix of eq.1 coefficients' draws
# beta_2(R X 3)  matrix of eq.2 coefficients' draws
# lambda (M_ind X 2 X R) array of mediator indicator coefficients' draws.
# tau (Y_ind X 2 X R) array of dependent variable indicator coefficients' draws.
# Each slice is one draw, where rows represent the indicator equation and columns are the coefficients
# All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively
# ssq_m_star(R X M_ind) Matrix of mediator indicator equations' coefficients' error variance draws
# ssq_y_star(R X Y_ind) Matrix of dependent variable indicator equation's coefficients' error variance draws
# mu_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)
# var_draw vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor)

MeasurementCont = function(Data, Prior, R, burnin){
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
  M = Data$M
  Y = Data$Y
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

  # R2 = R - burnin #account for burnin in the rest of computations

  #Running a Gibbs sampler regression of latent Y on latent M and X to store the MCMC draws for Bayes factor computation
  out = RuniregGibbsMulti(Data = list(y=post$Y, M = post$M, X=as.matrix(X)), Prior = list(ssq=1,A=diag(1/A_Y^2)), Mcmc = list(R=R-burnin,sigmasq=post$ssq_Y))

  return(list(beta_1 = cbind(post$beta_0_M,post$beta_1), beta_2 = cbind(post$beta_0_Y,post$beta_2,post$beta_3),
              ssq_M=post$ssq_M, ssq_Y = post$ssq_Y, Mdraw = post$M, Ydraw = post$Y,
              ssq_m_star = post$ssq_m_star, ssq_y_star = post$ssq_y_star,
              lambda = post$lambda, tau = post$tau,
              mu_draw = out$mubeta[,3], var_draw = out$IR[3,3,]))
  # return(list(beta_1 = cbind(post$beta_0_M,post$beta_1[1:R2]), beta_2 = cbind(post$beta_0_Y,post$beta_2,post$beta_3)[1:R2,],
  #             ssq_M=post$ssq_M[1:R2], ssq_y = post$ssq_Y[1:R2], Mdraw = post$M, Ydraw = post$Y,
  #             ssq_m_star = post$ssq_m_star[1:R2,], ssq_y_star = post$ssq_y_star[1:R2,],
  #             lambda = post$lambda[1:R2,,], tau = post$tau[1:R2,,],
  #             mu_draw = out$mubeta[,3], var_draw = out$IR[3,3,]))
}
