rbprobitGibbs_me=
function(Data,Prior,Mcmc)
{
#
# purpose:
#   draw from posterior for binary probit using Gibbs Sampler
#
# Arguments:
#   Data - list of X,y
#     X is nobs x nvar, y is nobs vector of 0,1
#   Prior - list of A, betabar
#     A is nvar x nvar prior preci matrix
#     betabar is nvar x 1 prior mean
#   Mcmc
#     R is number of draws
#     keep is thinning parameter
#     nprint - print estimated time remaining on every nprint'th draw
#
# Output:
#   list of betadraws
#
# Model:   y = 1 if  w=Xbeta + e   > 0  e ~N(0,1)
#
# Prior:   beta ~ N(betabar,A^-1)
#
#
#
# check arguments
#
if(missing(Data)) {stop("Requires Data argument -- list of y and X")}
    if(is.null(Data$X)) {stop("Requires Data element X")}
    X=Data$X
    if(is.null(Data$y)) {stop("Requires Data element y")}
    y=Data$y
nvar=ncol(X)
nobs=length(y)
#
# check data for validity
#
if(length(y) != nrow(X) ) {stop("y and X not of same row dim")}
if(sum(unique(y) %in% c(0:1)) < length(unique(y))) {stop("Invalid y, must be 0,1")}
#
# check for Prior
#
if(missing(Prior))
   { betabar=c(rep(0,nvar)); A=.01*diag(nvar)}
else
   {
    if(is.null(Prior$betabar)) {betabar=c(rep(0,nvar))}
       else {betabar=Prior$betabar}
    if(is.null(Prior$A)) {A=.01*diag(nvar)}
       else {A=Prior$A}
   }
#
# check dimensions of Priors
#
if(ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar)
   {stop(paste("bad dimensions for A",dim(A)))}
if(length(betabar) != nvar)
   {stop(paste("betabar wrong length, length= ",length(betabar)))}
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
    }
#
# # print out problem
# #
# cat(" ", fill=TRUE)
# cat("Starting Gibbs Sampler for Binary Probit Model",fill=TRUE)
# cat("   with ",length(y)," observations",fill=TRUE)
# cat("Table of y Values",fill=TRUE)
# print(table(y))
# cat(" ", fill=TRUE)
# cat("Prior Parms: ",fill=TRUE)
# cat("betabar",fill=TRUE)
# print(betabar)
# cat("A",fill=TRUE)
# print(A)
# cat(" ", fill=TRUE)
# cat("MCMC parms: ",fill=TRUE)
# cat("R= ",R," keep= ",keep," nprint= ",nprint,fill=TRUE)
# cat(" ",fill=TRUE)

beta=c(rep(0,nvar))
sigma=c(rep(1,nrow(X)))
root=chol(chol2inv(chol((crossprod(X,X)+A))))
Abetabar=crossprod(A,betabar)
        #a=ifelse(y == 0,-100, 0)
        #b=ifelse(y == 0, 0, 100)
        above=ifelse(y == 0, 1, 0)
        trunpt=0

###################################################################
# Keunwoo Kim
# 08/05/2014
# Edited by A.Laghaie 2020
###################################################################
draws=rbprobitGibbs_rcpp_me(y,X,Abetabar,root,beta,sigma,trunpt,above,R,keep,nprint)
###################################################################

return(draws)
}
