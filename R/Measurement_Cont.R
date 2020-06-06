#' Measurement_Cont
#'
#' @param Data asd
#' @param Prior asd
#' @param Mcmc asd
#'
#' @return asd
#' @export
#'
Measurement_Cont = function(Data, BF){
  if(missing(BF))
  { R = 5000; burnin = 3000; A_M = rep(10,2); A_Y = c(10,10,1);}
  else
  {
    if(is.null(BF$R)) {R = 5000}
    else {R = BF$R}
    if(is.null(BF$burnin)) {burnin = 3000}
    else {burnin = BF$burnin}
    if(is.null(BF$A_M)) {A_M = rep(10,2)}
    else {A_M = sqrt(BF$A_M)}
    if(is.null(BF$A_Y)) {A_Y = c(10,10,1)}
    else {A_Y = sqrt(BF$A_Y)}
  }

  X = Data$X
  M = Data$M
  Y = Data$Y
  N = length(Data$X);
  evidence = NULL


  m_star = Data$m_star
  y_star = as.matrix(Data$y_star)
  m_ind =  dim(m_star)[2];
  y_ind =  dim(y_star)[2];


stanfit = rstan::sampling(object = stanmodels$Measurement_Multi,
               data = list(n=N,M_ind = m_ind, Y_ind = y_ind,A_M=A_M, A_Y=A_Y, X = X, m_star=t(m_star),y_star=t(y_star)),
               pars = c("beta_1","beta_2","beta_3","ssq_M","ssq_Y","ssq_m_star","ssq_y_star","beta_0_Y","beta_0_M","lambda","tau","M","Y"),
               chains = 1,iter = R, warmup = burnin)

out.multi = rstan::extract(stanfit)

return(out.multi)
}
