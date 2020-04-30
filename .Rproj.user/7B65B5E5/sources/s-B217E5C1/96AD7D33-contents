#' BFSD
#'
#' @param Data Data
#' @param post Posterior
#' @param prior Prior
#' @param model Model
#' @param burnin Burnin
#' @param M MMMMM
#'
#' @return  HAHAHAHA
#' @export
#'
BFSD = function(Data,post,prior,model,burnin,M){

  #computes the Bayes Factor using Savage Dicky approxmiation
  if(missing(prior)){
    A = .01;          #A is the precision of beta_3
  }
  else
  {
    if(is.null(prior$A)) {A = .01}
    else {A = prior$A}
  }

  lik=unnorm_post=num=0;
  nobs = dim(Data$X)[1]
  if(is.vector(Data$X)) nobs=length(Data$X)

  ######################################
  ######################################


  ## For simple mediation and models with continuous indicators
  if(model=='cont'){  #(model=='measurement_SD'){

    R =  length(post$ssq_M)
    k=dim(as.matrix(post$beta_1))[2] #for stan computations k is =1 in unidimensional X case, and >1 for multidimensional x, for runiregGibbs (since we estimate intercept, k will be =2 and higher)
    Data$X = cbind(rep(1,nobs),Data$X)
    outer= rep(0,(R-burnin)) #the outer sum

    #computing the marginal posterior of beta_3
    for(r in (burnin+1):R){
        outer[r-burnin] = dnorm(x = 0, mean = post$mu_beta_2[r,k+1], sd = sqrt(post$IR_beta_2[k+1,k+1,r]), log = T)
    }

    #Numerically Stably computing the log of the outer sum
    outermax = max(outer)
    outer = outer[outer != max(outer)] #removes the maximum element from the vector
    lik = (outermax + log(1 + sum(exp(outer - outermax))))  #the outer sum in the log-scale
    lik = -log((R-burnin)) + lik

    lik = lik - dnorm(x = 0, mean = 0, sd = sqrt(1/A), log = T)
    lik = -lik #only to adjust to the function output
  }

  # if(model=='probit_SD'){
  #
  #   R =  dim(post$betadraw)[1]
  #   k=dim(as.matrix(post$betadraw))[2]
  #   Data$X = cbind(rep(1,nobs),Data$X)
  #   outer= rep(0,(R-burnin)) #the outer sum
  #
  #   #computing the marginal posterior of beta_3
  #   for(r in (burnin+1):R)  outer[r-burnin] = dnorm(x = 0, mean = post$mu_beta[r,k],
  #                                                   sd = sqrt(post$IR[k,k,r]), log = T)
  #
  #   #Numerically Stably computing the log of the outer sum
  #   outermax = max(outer)
  #   outer = outer[outer != max(outer)] #removes the maximum element from the vector
  #   lik = (outermax + log(1 + sum(exp(outer - outermax))))  #the outer sum in the log-scale
  #   lik = -log((R-burnin)) + lik
  #
  #   lik = lik - dnorm(x = 0, mean = 0, sd = sqrt(1/A), log = T)
  #   lik = -lik #only to adjust to the function output
  # }

  ## For models with categorical indicators
  if(model=='cat'){  #(model=='ordered'){

    R =  dim(post$betadraw)[1]
    k= 3; #dim(post$beta_2_draw)[2]
    outer= rep(0,(R-burnin)) #the outer sum

    #computing the marginal posterior of beta_3
    for(r in (burnin+1):R)  outer[r-burnin] = dnorm(x = 0, mean = post$mubeta_2_draw[r,k],
                                                    sd = sqrt(post$varbeta_2_draw[k,k,r]), log = T)

    #Numerically Stably computing the log of the outer sum
    outermax = max(outer)
    outer = outer[outer != max(outer)] #removes the maximum element from the vector
    lik = (outermax + log(1 + sum(exp(outer - outermax))))  #the outer sum in the log-scale
    lik = -log((R-burnin)) + lik

    lik = lik - dnorm(x = 0, mean = 0, sd = sqrt(1/A), log = T)
    lik = -lik #only to adjust to the function output
  }

  # if(model=='cat_Y'){
  #
  #   R =  dim(post$betadraw)[1]
  #   k = dim(post$betadraw)[2]
  #   outer= rep(0,(R-burnin)) #the outer sum
  #
  #   #computing the marginal posterior of beta_3
  #   for(r in (burnin+1):R)  outer[r-burnin] = dnorm(x = 0, mean = post$mubetadraw[r,k],
  #                                                   sd = sqrt(post$varbetadraw[k,k,r]), log = T)
  #
  #   #Numerically Stably computing the log of the outer sum
  #   outermax = max(outer)
  #   outer = outer[outer != max(outer)] #removes the maximum element from the vector
  #   lik = (outermax + log(1 + sum(exp(outer - outermax))))  #the outer sum in the log-scale
  #   lik = -log((R-burnin)) + lik
  #
  #   lik = lik - dnorm(x = 0, mean = 0, sd = sqrt(1/A), log = T)
  #   lik = -lik #only to adjust to the function output
  # }

  return(-lik)
}
