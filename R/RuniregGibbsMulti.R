RuniregGibbsMulti=
  function(Data,Prior,Mcmc)
  {
    # Rcpp::sourceCpp('runiregGibbs_rcpp_multi.cpp')
    #
    # revision history:
    #          P. Rossi 1/17/05
    #          3/07 added classes
    #          W. Taylor 4/15 - added nprint option to MCMC argument
    # Purpose:
    #   perform Gibbs iterations for Univ Regression Model using
    #     prior with beta, sigma-sq indep
    #
    # Arguments:
    #   Data -- list of data
    #           y,X
    #   Prior -- list of prior hyperparameters
    #     betabar,A      prior mean, prior precision
    #     nu, ssq        prior on sigmasq
    #   Mcmc -- list of MCMC parms
    #     sigmasq=initial value for sigmasq
    #     R number of draws
    #     keep -- thinning parameter
    #     nprint - print estimated time remaining on every nprint'th draw
    #
    # Output:
    #   list of beta, sigmasq
    #
    # Model:
    #   y = Xbeta + e  e ~N(0,sigmasq)
    #          y is n x 1
    #          X is n x k
    #          beta is k x 1 vector of coefficients
    #
    # Priors:  beta ~ N(betabar,A^-1)
    #          sigmasq ~ (nu*ssq)/chisq_nu
    #
    #
    # check arguments
    #
    if(missing(Data)) {stop("Requires Data argument -- list of y and X")}
    if(is.null(Data$X)) {stop("Requires Data element X")}
    X=Data$X
    if(is.null(Data$M)) {stop("Requires Data element y")}
    M=Data$M
    if(is.null(Data$y)) {stop("Requires Data element y")}
    y=Data$y
    nvar= 3  # our X has 3 columns: intercept, X and M
    nobs=ncol(y)
    #
    # check data for validity
    #
    if(nobs != nrow(X) ) {stop("length(y) ne nrow(X)")}
    #
    # check MCMC argument
    #
    if(missing(Mcmc)) {stop("requires Mcmc argument")}
    else
    {
      if(is.null(Mcmc$R))
      {stop("requires Mcmc element R")} else {R=Mcmc$R}
      if(is.null(Mcmc$keep)) {keep=1} else {keep=Mcmc$keep}
      if(is.null(Mcmc$nprint)) {nprint=100} else {nprint=Mcmc$nprint}
      if(nprint<0) {stop('nprint must be an integer greater than or equal to 0')}
      if(is.null(Mcmc$sigmasq)) {sigmasq=stats::var(y)} else {sigmasq=Mcmc$sigmasq}
    }
    #
    # check for Prior
    #
    if(missing(Prior))
    { betabar=c(rep(0,nvar)); A=.01*diag(nvar); nu=3; ssq=stats::var(y); betafix=F; sigmafix=F; betavalue=matrix(double(R*nvar),ncol=nvar); sigmavalue=rep(0,R)}
    else
    {
      if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))}
      else {betabar=Prior$betabar}
      if(is.null(Prior$A)) {A=.01*diag(nvar)}
      else {A=Prior$A}
      if(is.null(Prior$nu)) {nu=3}
      else {nu=Prior$nu}
      if(is.null(Prior$ssq)) {ssq=stats::var(y)}
      else {ssq=Prior$ssq}
      if(is.null(Prior$betafix)) {betafix=F}   #betafix is the indicator for setting the intercept to 0 and the beta coefficient to 1 (for estimating M in the measurement error equation )
      else {betafix=Prior$betafix}
      if(is.null(Prior$sigmafix)) {sigmafix=F}   #betafix is the indicator for setting the intercept to 0 and the beta coefficient to 1 (for estimating M in the measurement error equation )
      else {sigmafix=Prior$sigmafix}
      if(is.null(Prior$betavalue)) {betavalue=matrix(double(R*nvar),ncol=nvar)}   #betavalue is the fixed value of beta (in case betafix=T )
      else {betavalue=Prior$betavalue}
      if(is.null(Prior$sigmavalue)) {sigmavalue=rep(0,2)}   #sigmavalue is the fixed value of sigma (in case sigmafix=T )
      else {sigmavalue=Prior$sigmavalue}
    }
    #
    # check dimensions of Priors
    #
    if(ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar)
    {stop(paste("bad dimensions for A",dim(A)))}
    if(length(betabar) != nvar)
    {stop(paste("betabar wrong length, length= ",length(betabar)))}
    #
    # print out problem
    #
    # cat(" ", fill=TRUE)
    # cat("Starting Gibbs Sampler for Univariate Regression Model",fill=TRUE)
    # cat("  with ",nobs," observations",fill=TRUE)
    # cat(" ", fill=TRUE)
    # cat("Prior Parms: ",fill=TRUE)
    # cat("betabar",fill=TRUE)
    # print(betabar)
    # cat("A",fill=TRUE)
    # print(A)
    # cat("nu = ",nu," ssq= ",ssq,fill=TRUE)
    # cat(" ", fill=TRUE)
    # cat("MCMC parms: ",fill=TRUE)
    # cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
    # cat(" ",fill=TRUE)

    ###################################################################
    # Keunwoo Kim 08/05/2014
    #edited by A. Laghaie 2018 for the estimation of the measurement error
    ###################################################################
    draws = RuniregGibbsMultiCpp(y,M, t(X), betabar, A, nu, ssq, sigmasq, R, keep, nprint, betafix, sigmafix, betavalue, sigmavalue)
    ###################################################################

    # attributes(draws$betadraw)$class=c("bayesm.mat","mcmc")
    # attributes(draws$betadraw)$mcpar=c(1,R,keep)
    # attributes(draws$sigmasqdraw)$class=c("bayesm.mat","mcmc")
    # attributes(draws$sigmasqdraw)$mcpar=c(1,R,keep)
    # #Moments for calculating Bayes factor
    # attributes(draws$mubetadraw)$class=c("bayesm.mat","mcmc")
    # attributes(draws$mubetadraw)$mcpar=c(1,R,keep)
    # attributes(draws$IRdraw)$class=c("bayesm.mat","mcmc")
    # attributes(draws$IRdraw)$mcpar=c(1,R,keep)
    # attributes(draws$nudraw)$class=c("bayesm.mat","mcmc")
    # attributes(draws$nudraw)$mcpar=c(1,R,keep)
    # attributes(draws$ssqdraw)$class=c("bayesm.mat","mcmc")
    # attributes(draws$ssqdraw)$mcpar=c(1,R,keep)


    return(draws)
  }
