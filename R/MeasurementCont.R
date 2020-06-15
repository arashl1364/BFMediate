#' MeasurementCont
#'
#' @description No description yet
#'
#' @usage \code{MeasurementCont(Data, Prior, R, burnin)}
#'
#' @param Data list(X, m_tilde, y_tilde)
#' @param Prior list(A_M,A_Y)
#' @param R number of MCMC iterations, default = 10000
#'
#' @details
#' \subsection{Model}
#' \tabular{ll}{
#' M = beta_0M + Xbeta_1 + U_M  \tab \emph{[eq.1]} \cr
#' Y = beta_0Y + Mbeta_2 + Xbeta_3 + U_Y \tab \emph{[eq.2]} \cr
#' }
#' Indicator equations:
#' \tabular{lcl}{
#' m*_1    \tab = \tab M + U_m*_1 \cr
#' ˜m_1   \tab = \tab  OrdProbit(m*_1,C_m_1) \cr
#'  m*_2     \tab = \tab lambda_01 + M + U_m*_2 \cr
#'  ˜m_2  \tab = \tab OrdProbit(m*_2,C_m_2) \cr
#'  ... \tab  \tab  \cr
#' m*_k   \tab =  \tab lambda_0k-1 + M + U_m*_k \cr
#' ˜m_k   \tab = \tab OrdProbit(m*_k,C_m_k) \cr
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
#' \subsection{Argument Details}
#' \code{Data = list(X, m_tilde, y_tilde)}
#' \tabular{ll}{
#' \code{X(N x 1) } \tab treatment variable vector \cr
#' \code{m_tilde(N x M_ind) } \tab ediator indicators' matrix  \cr
#' \code{y_tilde(N x Y_ind) } \tab dependent variable indicators' matrix \cr
#' }
#'
#' \code{Prior = list(A_M,A_Y)} \emph{[optional]}
#' \tabular{ll}{
#' \code{A_M }   \tab vector of coefficients' prior variances of eq.1, default = rep(100,2) \cr
#' \code{A_Y }   \tab vector of coefficients' prior variances of eq.2, default = c(100,100,1) \cr
#' }
#'
#' @return
#' \tabular{ll}{
#' \code{beta_1(R X 2) } \tab  matrix of eq.1 coefficients' draws \cr
#' \code{beta_2(R X 3) } \tab  matrix of eq.2 coefficients' draws \cr
#' \code{lambda (M_ind X 2 X R) } \tab array of mediator indicator coefficients' draws. Each slice is one draw, where rows represent the indicator equation and columns are the coefficients. All Slope coefficients as well as intercept of the first equation are fixed to 1 and 0 respectively. \cr
#' \code{ssq_m_star(R X M_ind)} \tab  Matrix of mediator indicator equations' coefficients' error variance draws \cr
#' \code{ssq_y_star(R X Y_ind) } \tab  Matrix of dependent variable indicator equations' coefficients' error variance draws \cr
#' \code{mu_draw } \tab  vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor) \cr
#' \code{var_draw } \tab  vector of means of MCMC draws of the direct effect (used in BFSD to compute Bayes factor) \cr
#' }
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

  return(list(beta_1 = cbind(post$beta_0_M,post$beta_1), beta_2 = cbind(post$beta_0_Y,post$beta_2,post$beta_3),
              ssq_M=post$ssq_M, ssq_Y = post$ssq_Y, Mdraw = post$M, Ydraw = post$Y,
              ssq_m_star = post$ssq_m_star, ssq_y_star = post$ssq_y_star,
              lambda = post$lambda, tau = post$tau,
              mu_draw = out$mubeta[,3], var_draw = out$IR[3,3,]))
}
