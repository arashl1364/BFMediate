#' Bayes Factor Power Analysis
#'
#' @usage BF_power_analysis(pars, Model)
#'
#' @param pars  list(N, beta_M, beta_Y, sigma_M, sigma_Y, R, burnin, reps) for "Simple", list(N, beta_M, beta_Y, ssq_m_star, lambda, Sigma_Y, M_ind, Mcut, cutoff_M, R, burnin, reps) for "MCat", list(N, beta_M, beta_Y, ssq_y_star, tau, sigma_M, Y_ind, Ycut, cutoff_Y, R, burnin, reps)  for "YCat", and list(N, beta_M, beta_Y, ssq_m_star, ssq_y_star, lambda, Sigma_Y, M_ind, Mcut, cutoff_M, R, burnin, reps)  for "MYCat"
#' @param Model can be either "Simple", "MCat", "YCat", "MYCat". 
#'
#' @details
#' ## Model
#' For further information about the models, see
#'
#' * \link[BFMediate]{PartialMed} for "Simple"
#' * \link[BFMediate]{MeasurementMCat} for "MCat"
#' * \link[BFMediate]{MeasurementYCat} for "YCat"
#' * \link[BFMediate]{MeasurementMYCat} for "MYCat"
#'
#' @return
#' ## The function returns the Bayes factor distribution obtained from "reps" datasets generated with the input parameter values
#'
#' @export
#' @examples
#'pars = NULL
#'pars$M_ind = 2
#'pars$Y_ind = 2
#'pars$Mcut = pars$Ycut = 8
#'pars$N = 50
#'pars$beta_M = c(.5,1)
#'pars$beta_Y = c(.7, 1.5, 0)
#'pars$ssq_m_star = c(.5,.3)
#'pars$lambda = c(0,-.5) #the intercepts for the latent M indicators w. measurement error
#'pars$ssq_y_star = c(.2,.2)
#'pars$tau = c(0,-.5)
#'pars$cutoff_M = matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
#'                       -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= pars$Mcut, byrow = TRUE)
#'pars$cutoff_Y =  matrix(c(-100, 0, 1.6, 2, 2.2, 3.3, 6,  100,
#'                          -100, 0, 1, 2, 3, 4, 5, 100) ,ncol= pars$Ycut, byrow = TRUE)
#'pars$R = 10000
#'pars$burnin = 2000
#'pars$reps = 100  
#'
#'BF_Power_analysis(pars = pars, model = "MYCat")    # The analysis takes several minutes
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

BF_Power_analysis = function(pars,model){

simPartialMed = function(beta_M,beta_Y, sigma_M, sigma_Y,N,X) {
  eps_M = rnorm(N)*sigma_M # generate errors for M (independent)
  eps_Y = rnorm(N)*sigma_Y # generate errors for Y (independent)
  M = beta_M[1] + beta_M[2] * X + eps_M # generate latent mediator M
  Y = beta_Y[1] + beta_Y[2] * M + beta_Y[3] * X + eps_Y # generate dependent variable
  list(X = X, M = M, Y = Y)
}
SimMeasurementMCat = function(X, beta_M, cutoff_M, beta_Y, Sigma_Y, M_ind, lambda, ssq_m_star){
  
  nobs = dim(X)[1]
  m_star = m_tilde = matrix(double(nobs*M_ind), ncol = M_ind)
  
  M = beta_M[1] + beta_M[2] * X + rnorm(nobs)  #cbind(rep(1,nobs),X)%*%beta_M + rnorm(nobs)
  
  for(i in 1: M_ind){
    m_star[,i] = lambda[i] + M + sqrt(ssq_m_star[i])*rnorm(nobs);
    m_tilde[,i] = cut(m_star[,i], br = cutoff_M[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
  }
  
  Y = beta_Y[1] + beta_Y[2] * M + beta_Y[3] * X + rnorm(nobs)*Sigma_Y
  #cbind(rep(1,nobs),cbind(M,X))%*%beta_Y + rnorm(nobs)
  
  return(list(Y = Y, M = M, m_tilde = m_tilde, X = X,
              beta_M = beta_M, beta_Y = beta_Y,
              lambda = lambda, ssq_m_star = ssq_m_star, m_star = m_star, cutoff_M = cutoff_M,
              k_M=dim(cutoff_M)[2]-1, M_ind=M_ind))
}
SimMeasurementYCat = function(X, beta_M, beta_Y, sigma_M, cutoff_Y, Y_ind, tau, ssq_y_star){
  
  nobs = dim(X)[1]
  y_tilde = y_star = matrix(double(nobs*Y_ind), ncol = Y_ind)
  
  M = beta_M[1] + beta_M[2] * X + rnorm(nobs) * sigma_M
  Y = beta_Y[1] + beta_Y[2] * M + beta_Y[3] * X + rnorm(nobs)
  
  for(i in 1: Y_ind){
    y_star[,i] = tau[i] + Y + sqrt(ssq_y_star[i])*rnorm(nobs);
    y_tilde[,i] = cut(y_star[,i], br = cutoff_Y[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
  }
  
  return(list(Y = Y, M = M, y_tilde = y_tilde, X = X,
              beta_M = beta_M,
              k_Y=dim(cutoff_Y)[2]-1, beta_Y = beta_Y, tau = tau,
              ssq_y_star = ssq_y_star, y_star = y_star, cutoff_Y = cutoff_Y,
              Y_ind=Y_ind))
}
SimMeasurementMYCat = function(X, beta_M, cutoff_M, beta_Y, cutoff_Y, M_ind, Y_ind,
                               lambda, tau, ssq_m_star, ssq_y_star){
  
  nobs = dim(X)[1]
  m_tilde = m_star = matrix(double(nobs*M_ind), ncol = M_ind)
  y_tilde = y_star = matrix(double(nobs*Y_ind), ncol = Y_ind)
  
  M = cbind(rep(1,nobs),X)%*%beta_M + rnorm(nobs)
  for(i in 1: M_ind){
    m_star[,i] = lambda[i] + M + sqrt(ssq_m_star[i])*rnorm(nobs);
    m_tilde[,i] = cut(m_star[,i], br = cutoff_M[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
  }
  
  Y = cbind(rep(1,nobs),cbind(M,X))%*%beta_Y + rnorm(nobs)
  for(i in 1: Y_ind){
    y_star[,i] = tau[i] + Y + sqrt(ssq_y_star[i])*rnorm(nobs);
    y_tilde[,i] = cut(y_star[,i], br = cutoff_Y[i,], right=TRUE, include.lowest = TRUE, labels = FALSE)
  }
  
  return(list(Y = Y, M = M, y_tilde = y_tilde, m_tilde = m_tilde, X = X,
              k_M=dim(cutoff_M)[2]-1, beta_M = beta_M, lambda = lambda,
              ssq_m_star = ssq_m_star, m_star = m_star, cutoff_M = cutoff_M,
              k_Y=dim(cutoff_Y)[2]-1, beta_Y = beta_Y, tau = tau, ssq_y_star = ssq_y_star,
              y_star = y_star, cutoff_Y = cutoff_Y, M_ind=M_ind, Y_ind=Y_ind))
}

M_ind = pars$M_ind
Y_ind = pars$Y_ind
Mcut = pars$Mcut
Ycut = pars$Ycut
beta_M = pars$beta_M
beta_Y = pars$beta_Y
sigma_M = pars$sigma_M
sigma_Y = pars$sigma_Y
Sigma_Y = pars$Sigma_Y
ssq_m_star = pars$ssq_m_star
lambda = pars$lambda  #the intercepts for the latent M indicators w. measurement error
ssq_y_star = pars$ssq_y_star
tau = pars$tau
cutoff_M = pars$cutoff_M
cutoff_Y = pars$cutoff_Y
N = pars$N   # number of observations    
R = pars$R; burnin = pars$burnin 
reps = pars$reps

# k=length(beta_1)-1
reps = reps
A_M=rep(100,2); 
A_Y=c(100,100,1)


# CI_single = CI_multi = BF_single = BF_multi = matrix(double(reps*length(N)),ncol=length(N)) # columns 1,2, and 3 of BF contains BF values of replications of sample sizes 50,200, and 2000 respectively
# beta_3_post = matrix(double((R-burnin)*length(N)),ncol=length(N))
BF = rep(0,reps)


for(i in 1:reps){
  
  if(model=="Simple"){
    ##### Single #####   
    # Generate data based on parameters
    X = as.matrix(rbinom(N,1,c(.5,.5))) 
    Data = simPartialMed(beta_M,beta_Y,sigma_M,sigma_Y,N,X)
    invisible(capture.output(out <- PartialMed(Data=Data, Prior = list(A_M=A_M, A_Y=A_Y), R = R)))
    BF[i] = exp(BFSD(Post = out , Prior = A_Y[3], burnin = burnin))
    
  }
  
  if(model=="MCat"){
    ##### multiple indicator with M and Y categorized #####   
    # Generate data based on parameters
    X = as.matrix(rbinom(N,1,c(.5,.5))) 
    DataMCat = SimMeasurementMCat(X, beta_M, cutoff_M, beta_Y, Sigma_Y, M_ind, lambda, ssq_m_star)
    Mcut = max(DataMCat$m_tilde) + 1
    Data = list(X=cbind(rep(1,length(DataMCat$X)),DataMCat$X), m_tilde=as.matrix(DataMCat$m_tilde),
                Y= as.matrix(DataMCat$Y) ,k=Mcut-1, M_ind=dim(DataMCat$m_tilde)[2])
    Prior = list(A_M = A_M, A_Y = A_Y)
    invisible(capture.output(out <- MeasurementMCat(Data=Data, Prior=Prior, R=R)))
    BF[i] = exp(BFSD(Post = out , Prior = A_Y[3], burnin = burnin))
    
  }
  
  if(model=="YCat"){
    ##### multiple indicator with M and Y categorized #####   
    # Generate data based on parameters
    X = as.matrix(rbinom(N,1,c(.5,.5))) 
    DataYCat = SimMeasurementYCat(X, beta_M, beta_Y, sigma_M, cutoff_Y, Y_ind, tau, ssq_y_star)
    Ycut = max(as.matrix(DataYCat$y_tilde)[,1]) +1
    Data = list(X=cbind(rep(1,length(DataYCat$X)),DataYCat$M,DataYCat$X), y = as.matrix(DataYCat$y_tilde),
                k=Ycut-1, Y_ind=dim(as.matrix(DataYCat$y_tilde))[2])
    Prior = list(A_M = A_M, A_Y = A_Y)
    invisible(capture.output(out <- MeasurementYCat(Data=Data, Prior=Prior, R=R)))
    BF[i] = exp(BFSD(Post = out , Prior = A_Y[3], burnin = burnin))
    
  }
  
  if(model=="MYCat"){
    ##### multiple indicator with M and Y categorized #####   
    # Generate data based on parameters
    X = as.matrix(rbinom(N,1,c(.5,.5))) 
    DataMYCat = SimMeasurementMYCat(X, beta_M, cutoff_M, beta_Y, cutoff_Y, M_ind,
                                    Y_ind, lambda, tau, ssq_m_star, ssq_y_star)
    Mcut = max(DataMYCat$m_tilde) + 1
    Ycut = max(DataMYCat$y_tilde) + 1
    Data = list(X=cbind(rep(1,length(DataMYCat$X)),DataMYCat$X), m_tilde=as.matrix(DataMYCat$m_tilde),
                y_tilde=as.matrix(DataMYCat$y_tilde), k_M = Mcut-1, k_Y=Ycut-1,
                M_ind=dim(as.matrix(DataMYCat$m_tilde))[2], Y_ind=dim(as.matrix(DataMYCat$y_tilde))[2])
    Prior = list(A_M = A_M, A_Y = A_Y)
    invisible(capture.output(out <- MeasurementMYCat(Data=Data, Prior=Prior, R=R)))
    BF[i] = exp(BFSD(Post = out , Prior = A_Y[3], burnin = burnin))
    
  }

  if(i/10==floor(i/10)){cat("Replication:",i,"\n")}
  
}
## Plotting BF density
BF_quant = quantile(BF,probs = c(.025,.975))  
dense_BF = density(BF)
ref_len = max(dense_BF$y)/50 #length of the reference lines for mean and quantiles
BF_range = rep(0,2)
BF_range[1] = min(dense_BF$x) #- 1
BF_range[2] = max(dense_BF$x) #+ 1

df = data.frame(N=rep(pars$N,pars$N),BF = BF)
quartz();
ggplot2::ggplot(df, aes(x=BF)) + 
  # theme_bw() +
  geom_density(color = 'deeppink',fill='deeppink4') +
  geom_segment(aes(x = BF_quant[1] , y = 0 , xend = BF_quant[1] , yend=ref_len),color = "yellow")  +
  geom_segment(aes(x = BF_quant[2] , y = 0 , xend = BF_quant[2] , yend=ref_len),color = "yellow") +
  geom_segment(aes(x = mean(BF) , y = 0 , xend = mean(BF) , yend=ref_len),color = "orange",size=1.5) +
  xlim(BF_range)
  
}
