#' Title
#'
#' @param Data asdf
#' @param post asdf
#' @param pars asdf
#' @param R asdf
#'
#' @return asdf
#' @export
#'
GenMoments=function(Data, post, pars, R){
  #Moment generator (and augmentation of M) of a latent mediation effect when there is a measurement error
  #and direct effect. Get the stan posterior as input and gives M and moments as output for BF computation
  #Data: X(Nxk),y(Nx1),m_star(Nx1)
  #Parameters:  M, beta_1, beta_2, ssq_M, ssq_m_star, ssq_y
  ######## NEW : Since the new Stan code estimates the intercepts, no need to demean the data ########

  #Prior and parameter initial values
  X = as.matrix(Data$X)

  if(missing(pars)){
    A = diag(.01,ncol(X)+2)
  }
  else
  {
    if(is.null(pars$A)) {A = diag(.01,ncol(X))}
    else {A = diag(pars$A,ncol(X)+2)}
  }

  N=dim(X)[1]
  k=dim(X)[2]
  draws=length(post$ssq_M)


  ##Posterior draws
  lambda = post$lambda
  tau = post$tau
  Mdraw = ydraw = matrix(double(R*N),ncol = N)
  beta_0_M = post$beta_0_M
  beta_0_Y = post$beta_0_Y
  beta_1 = post$beta_1;
  beta_2 = post$beta_2;
  beta_3 = post$beta_3;
  ssq_M = post$ssq_M;
  ssq_Y = post$ssq_Y;
  ssq_m_star = post$ssq_m_star
  ssq_y_star = post$ssq_y_star
  M = post$M
  Y = post$Y

  ##Moment draws
  mu_beta_2_draw = matrix(double(R*(k+2)),ncol = k+2)
  IR_beta_2_draw =  array(double(((k+2)^2)*R), dim = c(k+2,k+2,R))
  #

  itime=proc.time()[3]
  cat("MCMC Iteration (est time to end - min) ",fill=TRUE)


  # for(i in 1:R){
  #
  #
  #   # #draw beta_2,beta_3 ssq_y | y,M,X
  #
  #   out<-runiregGibbs_me(Data = list(y=Y[i,],X=cbind(rep(1,N),M[i,],X)),Prior=list(ssq=1), Mcmc = list(R=1,sigmasq=ssq_Y[i]))
  #
  #   ##Moments
  #   mu_beta_2 = out$mubeta
  #   IR_beta_2 = out$IR
  #   nu_ssq_y = out$nu
  #   S_ssq_y = out$S
  #
  #   #storing moments
  #
  #   mu_beta_2_draw[i,] = mu_beta_2
  #   IR_beta_2_draw[,,i] = IR_beta_2
  #
  #
  #   #       print time to completion and draw # every 100th draw
  #   #
  #   if(i%%100 == 0)
  #   {ctime=proc.time()[3]
  #   timetoend=((ctime-itime)/i)*(R-i)
  #   cat(" ",i," (",round(timetoend/60,1),")",fill=TRUE)
  #   #fsh()
  #   }
  # }
  ctime = proc.time()[3]
  cat('  Total Time Elapsed: ',round((ctime-itime)/60,2),'\n')

  out<-RuniregGibbsMulti(Data = list(y=Y, M = M, X=X), Prior = list(ssq=1,A=A), Mcmc = list(R=R,sigmasq=ssq_Y))

  mu_beta_2_draw = out$mubeta
  IR_beta_2_draw = out$IR

  return(list(beta_1=cbind(beta_0_M,beta_1[1:R]),beta_2=cbind(beta_0_Y,beta_2,beta_3)[1:R,], ssq_M=ssq_M[1:R] ,ssq_y=ssq_Y[1:R], Mdraw=M, ydraw = Y,
              ssq_m_star=ssq_m_star[1:R,], ssq_y_star=ssq_y_star[1:R,], lambda= lambda[1:R,,], tau=tau[1:R,,],
              mu_beta_2=mu_beta_2_draw, IR_beta_2=IR_beta_2_draw))
}


