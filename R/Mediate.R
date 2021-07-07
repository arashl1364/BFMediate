#' Mediation Analysis and Bayes Factor Computation
#'
#' @usage Mediate(Data, Model, Prior, R, burnin)
#'
#' @param Data  list(X, M, Y) for "Simple", list(X, m_star, y_star) for "Cont", list(X, m_tilde, Y) for "MCat", list(X, M, y_tilde) for "YCat", and list(X, m_tilde, y_tilde) for "MYCat"
#' @param Model can be either "Simple", "Cont", "MCat", "YCat", "MYCat". In case of Simple, a simple partial mediation is estimated, Baron and Kenny (1986), and Preacher and Hayes (2004) proposed methods are also computed
#' @param Prior list(A_M,A_Y)
#' @param R number of MCMC iterations, default = 10000
#' @param burnin number of MCMC draws before the posterior is converged
#'
#' @details
#' ## Model
#' For Data arguments and Models, see
#'
#' * \link[BFMediate]{PartialMed} for "Simple"
#' * \link[BFMediate]{MeasurementMCat} for "MCat"
#' * \link[BFMediate]{MeasurementYCat} for "YCat"
#' * \link[BFMediate]{MeasurementMYCat} for "MYCat"
#'
#' \code{Prior = list(A_M,A_Y) }  *\[optional\]*
#'
#' \describe{
#'   \item{\code{A_M}}{vector of coefficients' prior variances of eq.1, default = rep(100,2)}
#'   \item{\code{A_Y}}{vector of coefficients' prior variances of eq.2, default = c(100,100,1)}
#' }
#'
#' @return
#' ## \code{BK = list(eq1, eq2, Indirect_se, FullMed)} (only for "Simple"!)
#' \describe{
#'   \item{eq1}{the summary of the eq.1 regression}
#'   \item{eq2}{the summary of the eq.2 regression}
#'   \item{Indirect_se}{the standard error of the indirect effect a la Sobel(1982)}
#'   \item{FullMed}{the significance test result for the direct effect}
#' }
#'
#' ## \code{PH = list(Indirect_mean, Indirect_CI, Direct_CI)}
#' \describe{
#'   \item{Indirect_mean}{the bootstrapped mean of the indirect effect}
#'   \item{Indirect_CI}{the bootstrapped 95% confidence interval of the indirect effect}
#'   \item{Direct_CI}{the bootstrapped 95% confidence interval of the direct effect}
#' }
#'
#' ## \code{list(evidence, Indirect_CI, Direct_CI, BF,...)} (For all the models)
#' \describe{
#'   \item{evidence}{the interpretation of the BF in terms of evidence in favor of full mediation according to Kass and Raftery (1995) }
#'   \item{Indirect_CI}{the Bayesian 95% HDI (confidence interval) of the indirect effect }
#'   \item{Direct_CI}{the Bayesian 95% HDI (confidence interval) of the direct effect}
#'   \item{BF}{the Bayes factor(BF_01) of the corresponding model (see Laghaie and Otter (2020))}
#' }
#'
#' For the rest of the values, see
#' * \link[BFMediate]{PartialMed} fpr "Simple"
#' * \link[BFMediate]{MeasurementCont} for "Cont"
#' * \link[BFMediate]{MeasurementMCat} for "MCat"
#' * \link[BFMediate]{MeasurementYCat} for "YCat"
#' * \link[BFMediate]{MeasurementMYCat} for "MYCat"
#' @export
#' @examples
#' simPartialMed = function(beta_M,beta_Y, sigma_M, sigma_Y,N,X) {
#' eps_M = rnorm(N)*sigma_M      # generate errors for M (independent)
#' eps_Y = rnorm(N)*sigma_Y      # generate errors for Y (independent)
#' M = beta_M[1] + beta_M[2] * X + eps_M # generate latent mediator M
#' Y = beta_Y[1] + beta_Y[2] * M + beta_Y[3] * X + eps_Y  # generate dependent variable
#' list(X = X, M = M, Y = Y)
#' }
#'
#' # Set up data generating parameters
#' N = 1000    # number of observations
#' sigma_M = .2^.5    # error std M
#' sigma_Y = .2^.5    # error std Y
#' beta_M = c(1, .3)   # beta_0M and beta_1
#' beta_Y = c(1, .5, .01)    # beta_0Y, beta_2, beta_3
#' X = rnorm(N,mean = 1,sd = 1)   # generate random X
#' # Generate data based on parameters
#' Data = simPartialMed(beta_M,beta_Y,sigma_M,sigma_Y,N,X)
#'
#' #Estimation
#' A_M = c(100,100);  # Prior variance for beta_0M, beta_1
#' A_Y = c(100,100,1) # Prior variance for beta_0Y, beta_2, beta_3
#' R = 2000
#' out = Mediate(Data = Data, Model = 'Simple',
#'                 Prior = list(A_M = A_M, A_Y = A_Y),R=5000, burnin = 3000)
#'
#' # Results
#' out$BK$FullMed
#' out$PH$Indirect_CI
#' colMeans(out$Simple$beta_2)
#' out$Simple$BF
### Description BFMediate estimates different mediation models and computes Bayes factors
# to test full mediation in them

### Arguments:
# Model can be either "Simple", "Cont", "MCat", "YCat", "MYCat".
# In case of Simple, a simple partial mediation is estimated,
# Baron and Kenny (1986), and Preacher and Hayes (2004) proposed methods are also computed
# Data  list(X, M, Y) for "Simple", list(X, m_star, y_star) for "Cont",
# list(X, m_tilde, Y) for "MCat", list(X, M, y_tilde) for "YCat",
# and list(X, m_tilde, y_tilde) for "MYCat"
# Prior list(A_M,A_Y)
# R number of MCMC iterations
# burnin number of MCMC draws before the posterior is converged

### Details:
## Model:
## For Data arguments and Models
## see PartialMed" ##KaiLong: Please make hyperlink to PartialMed, for others in this section do the same## , for "Simple"
## MeasurementCont for "Cont"
## MeasurementMCat for "MCat"
## MeasurementYCat for "YCat"
## MeasurementMYCat for "MYCat"
# Prior = list(A_M,A_Y) [optional]
# A_M vector of coefficients' prior variances of eq.1 (def: rep(100,2))
# A_Y vector of coefficients' prior variances of eq.2 (def: c(100,100,1))
# R number of MCMC iterations (def:10000)
# burnin number of MCMC iterations to be discarded from the draws

### Value:
# only for "Simple": BK = list(eq1, eq2, Indirect_se, FullMed)
# eq1 is the summary of the eq.1 regression
# eq2 is the summary of the eq.2 regression
# Indirect_se is the standard error of the indirect effect a la Sobel(1982)
# FullMed is the significance test result for the direct effect
# PH = list(Indirect_mean, Indirect_CI, Direct_CI)
# Indirect_mean is the bootstrapped mean of the indirect effect
# Indirect_CI is the bootstrapped 95% confidence interval of the indirect effect
# Direct_CI is the bootstrapped 95% confidence interval of the direct effect
# for all the models list(evidence, Indirect_CI, Direct_CI, BF,...)
# Indirect_CI is the Bayesian 95% HDI (confidence interval) of the indirect effect
# Direct_CI is the Bayesian 95% HDI (confidence interval) of the direct effect
# BF is the Bayes factor(BF_01) of the corresponding model (see Laghaie and Otter (2020) )
# evidence is the interpretation of the BF in terms of evidence in favor of full mediation according to Kass and Raftery (1995)
# for the rest of the values
## see PartialMed" ##KaiLong: Please make hyperlink to PartialMed, for others in this section do the same## , for "Simple"
## MeasurementCont for "Cont"
## MeasurementMCat for "MCat"
## MeasurementYCat for "YCat"
## MeasurementMYCat for "MYCat"
Mediate = function(Data, Model, Prior, R, burnin){  # BF){

  ############################################
  ## A. Laghaie 2019
  ############################################

  if(missing(Prior))
  { A_M = rep(100,2); A_Y = c(100,100,1);}
  else
  {
    if(is.null(Prior$A_M)) {A_M = rep(100,2)}
    else {A_M = Prior$A_M}
    if(is.null(Prior$A_Y)) {A_Y = c(100,100,1)}
    else {A_Y = Prior$A_Y}
  }
  if(missing(Model)) stop("Specify the Model of analysis")

  X = Data$X
  N = length(Data$X);
  evidence = NULL

  if(Model == 'Simple'){

    M = Data$M
    Y = Data$Y

    # B&K
    BK1 = summary(lm(formula = M ~ X))
    BK2 = summary(lm(formula = Y ~ M + X))
    BK.beta_1 = BK1$coefficients[2,1]
    BK.beta_1.sd = BK1$coefficients[2,2]
    BK.beta_2 = BK2$coefficients[2,1]
    BK.beta_2.sd = BK2$coefficients[2,2]
    BK.beta_3 = BK2$coefficients[3,1]
    BK.beta_3.sd = BK2$coefficients[3,2]
    BK.beta_3.pvalue = BK2$coefficients["X","Pr(>|t|)"]

    # B&K full mediation test
    std_res = ifelse(BK2$coefficients["X","Pr(>|t|)"]<.05,"Reject", "Fail to reject")

    #List of B&K results
    BK = list(eq1 = BK1, eq2 = BK2,
              # beta_1 = BK.beta_1, beta_2 = BK.beta_2, beta_3 = BK.beta_3,beta_1_se = BK.beta_1.sd, beta_2_se = BK.beta_2.sd, beta_3_se = BK.beta_3.sd,
              Indirect_se = sqrt(BK.beta_1^2*BK.beta_2.sd^2 + BK.beta_2^2*BK.beta_1.sd^2 + BK.beta_1.sd^2*BK.beta_2.sd^2),
              FullMed = std_res)


    # Preacher & Hayes bootstrapping
    boot = 5000
    b1 = b2 = b3 = rep(0,boot);
    for(i in 1:boot){
      samp = sample(length(X),replace = T)
      X_boot = X[samp]; M_boot = M[samp]; Y_boot = Y[samp]
      temp = summary(lm(formula = M_boot ~ X_boot))
      b1[i] = temp$coefficients[2,1]
      temp = summary(lm(formula = Y_boot ~ M_boot + X_boot))
      b2[i] = temp$coefficients[2,1]
      b3[i] = temp$coefficients[3,1]
    }
    PH.mean.Indirect = mean(b1*b2)
    PH.CI.Indirect = round(as.vector(quantile(b1*b2,probs=c(.025,.975))),2)
    PH.CI.Direct = round(as.vector(quantile(b3,probs=c(.025,.975))),2)

    PH = list(Indirect_mean = PH.mean.Indirect, Indirect_CI = PH.CI.Indirect, Direct_CI = PH.CI.Direct)


    # BF Simple
    Simple = PartialMed(Data = Data, R=R, Prior = list(A_M=A_M, A_Y=A_Y))
    beta_1 = Simple$beta_M[,2]
    beta_2 = Simple$beta_Y[,2]
    beta_3 = Simple$beta_Y[,3]
    BF.Simple = exp(BFSD(Post = Simple , Prior = A_Y[3], burnin = burnin))
    Bayes.CI.Indirect = round(as.vector(quantile(beta_1*beta_2,probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(beta_3,probs = c(.025,.975))),2)

    if(BF.Simple>1) evidence = ifelse(BF.Simple>100,"Decisive in favor of full mediation",
                                   ifelse(BF.Simple>10,"Strong in favor of full mediation",
                                          ifelse(BF.Simple>3.2,"Substantial in favor of full mediation","Not worth more than a bare mention")))

    if(BF.Simple<1) evidence = ifelse(1/BF.Simple>100,"Decisive against full mediation",
                                   ifelse(1/BF.Simple>10,"Strong against full mediation",
                                          ifelse(1/BF.Simple>3.2,"Substantial against full mediation","Not worth more than a bare mention")))

    # CI = as.character(ifelse((quantile(beta_3,probs = .025)>0) | (quantile(beta_3,probs = .975)<0),"Reject","Accept"))

    Simple$evidence = evidence
    Simple$Indirect_CI = Bayes.CI.Indirect
    Simple$Direct_CI = Bayes.CI.Direct
    Simple$BF = BF.Simple

    return(list(BK = BK, PH = PH, Simple = Simple))

    }

  ################################
  ##### Continous Data Multi #####
  ################################

  # BF multi
  if(Model == 'Cont'){

    m_star = Data$m_star
    y_star = as.matrix(Data$y_star)
    m_ind =  dim(m_star)[2];
    y_ind =  dim(y_star)[2];

    out = MeasurementCont(Data = Data, Prior = list(A_M = A_M, A_Y = A_Y),R=R, burnin = burnin)
    BF.LVM = exp(BFSD(Post = out , Prior = A_Y[3], burnin = 0)) #we already accounted for burnin in estimation
    Bayes.CI.Indirect = round(as.vector(quantile(out$beta_M[,2]*out$beta_Y[,2],probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(out$beta_Y[,3],probs = c(.025,.975))),2)

    if(BF.LVM>1) evidence = ifelse(BF.LVM>100,"Decisive in favor of full mediation",
                                   ifelse(BF.LVM>10,"Strong in favor of full mediation",
                                          ifelse(BF.LVM>3.2,"Substantial in favor of full mediation","Not worth more than a bare mention")))

    if(BF.LVM<1) evidence = ifelse(1/BF.LVM>100,"Decisive against full mediation",
                                   ifelse(1/BF.LVM>10,"Strong against full mediation",
                                          ifelse(1/BF.LVM>3.2,"Substantial against full mediation","Not worth more than a bare mention")))

    CI = as.character(ifelse((quantile(out$beta_Y[,3], probs = .025)>0) | (quantile(out$beta_Y[,3],probs = .975)<0),"Reject","Accept"))


    out$evidence = evidence
    out$Indirect_CI = Bayes.CI.Indirect
    out$Direct_CI = Bayes.CI.Direct
    out$BF = BF.LVM

    return(out)
    # return(list(Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
    #             out = out, BF.LVM = BF.LVM, CI = CI, evidence = evidence))

  }

  ################################
  #### Categorical Data Multi ####
  ################################

  # BF only M categorical
  if(Model == 'MCat'){

    Mcut = max(Data$m_tilde) +1
    Data_cat=list(X=cbind(rep(1,length(Data$X)),Data$X), m_tilde=as.matrix(Data$m_tilde), Y= as.matrix(Data$Y) ,k=Mcut-1, M_ind=dim(Data$m_tilde)[2])
    out = MeasurementMCat(Data=Data_cat, Prior = Prior, R=R) #rordprobitGibbs_me_M_multi_merr_cpp(Data=Data_cat, Mcmc=Mcmc)
    BF.LVM = exp(BFSD(Post = out , Prior = A_Y[3], burnin = burnin))
    Bayes.CI.Indirect = round(as.vector(quantile(out$beta_M[,2]*out$beta_Y[,2],probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(out$beta_Y[,3],probs = c(.025,.975))),2)

    if(BF.LVM>1) evidence = ifelse(BF.LVM>100,"Decisive in favor of full mediation",
                                   ifelse(BF.LVM>10,"Strong in favor of full mediation",
                                          ifelse(BF.LVM>3.2,"Substantial in favor of full mediation","Not worth more than a bare mention")))

    if(BF.LVM<1) evidence = ifelse(1/BF.LVM>100,"Decisive against full mediation",
                                   ifelse(1/BF.LVM>10,"Strong against full mediation",
                                          ifelse(1/BF.LVM>3.2,"Substantial against full mediation","Not worth more than a bare mention")))

    # CI = as.character(ifelse((quantile(out$beta_2[,3], probs = .025)>0) | (quantile(out$beta_2[,3],probs = .975)<0),"Reject","Accept"))

    out$evidence = evidence
    out$Indirect_CI = Bayes.CI.Indirect
    out$Direct_CI = Bayes.CI.Direct
    out$BF = BF.LVM

    return(out)
    # return(list(Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
    #             out = out, BF.LVM = BF.LVM, CI = CI, evidence = evidence))

  }
  # BF only Y categorical
  if(Model == 'YCat'){


    Ycut = max(as.matrix(Data$y_tilde)[,1]) +1
    Data_cat=list(X=cbind(rep(1,length(Data$X)),Data$M,Data$X), y = as.matrix(Data$y_tilde) ,k=Ycut-1, Y_ind=dim(as.matrix(Data$y_tilde))[2])
    Mcmc=list(R=R)
    out = MeasurementYCat(Data=Data_cat, Prior=Prior, R=R) #rordprobitGibbs_me_multi_merr_cpp(Data=Data_cat, Mcmc=Mcmc)
    BF.LVM = exp(BFSD(Post = out , Prior = A_Y[3], burnin = burnin))
    Bayes.CI.Indirect = round(as.vector(quantile(out$beta_M[,2]*out$beta_Y[,2], probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(out$beta_Y[,3],probs = c(.025,.975))),2)

    if(BF.LVM>1) evidence = ifelse(BF.LVM>100,"Decisive in favor of full mediation",
                                   ifelse(BF.LVM>10,"Strong in favor of full mediation",
                                          ifelse(BF.LVM>3.2,"Substantial in favor of full mediation","Not worth more than a bare mention")))

    if(BF.LVM<1) evidence = ifelse(1/BF.LVM>100,"Decisive against full mediation",
                                   ifelse(1/BF.LVM>10,"Strong against full mediation",
                                          ifelse(1/BF.LVM>3.2,"Substantial against full mediation","Not worth more than a bare mention")))

    # CI = as.character(ifelse((quantile(out$beta_2[,3], probs = .025)>0) | (quantile(out$beta_2[,3],probs = .975)<0),"Reject","Accept"))

    out$evidence = evidence
    out$Indirect_CI = Bayes.CI.Indirect
    out$Direct_CI = Bayes.CI.Direct
    out$BF = BF.LVM

    return(out)

    # return(list(Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
    #             out = out, BF.LVM = BF.LVM, CI = CI, evidence = evidence))

  }

  # BF both M and Y categorical
  if(Model == 'MYCat'){

    Mcut = max(Data$m_tilde) +1
    Ycut = max(Data$y_tilde) +1
    Data_cat=list(X=cbind(rep(1,length(Data$X)),Data$X), m_tilde=as.matrix(Data$m_tilde), y_tilde=as.matrix(Data$y_tilde), k_M = Mcut-1, k_Y=Ycut-1, M_ind=dim(as.matrix(Data$m_tilde))[2], Y_ind=dim(as.matrix(Data$y_tilde))[2])
    out = MeasurementMYCat(Data=Data_cat, Prior=Prior, R=R) #Mediation_Ordered_Multi_Merr(Data=Data_cat, Mcmc=Mcmc)
    BF.LVM = exp(BFSD(Post = out , Prior = A_Y[3], burnin = burnin))
    Bayes.CI.Indirect = round(as.vector(quantile(out$beta_M[,2]*out$beta_Y[,2],probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(out$beta_Y[,3],probs = c(.025,.975))),2)

    if(BF.LVM>1) evidence = ifelse(BF.LVM>100,"Decisive in favor of full mediation",
                                   ifelse(BF.LVM>10,"Strong in favor of full mediation",
                                          ifelse(BF.LVM>3.2,"Substantial in favor of full mediation","Not worth more than a bare mention")))

    if(BF.LVM<1) evidence = ifelse(1/BF.LVM>100,"Decisive against full mediation",
                                   ifelse(1/BF.LVM>10,"Strong against full mediation",
                                          ifelse(1/BF.LVM>3.2,"Substantial against full mediation","Not worth more than a bare mention")))

    # CI = as.character(ifelse((quantile(out$beta_2[,3], probs = .025)>0) | (quantile(out$beta_2[,3],probs = .975)<0),"Reject","Accept"))

    out$evidence = evidence
    out$Indirect_CI = Bayes.CI.Indirect
    out$Direct_CI = Bayes.CI.Direct
    out$BF = BF.LVM

    return(out)

    # return(list(Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
    #             out = out, BF.LVM = BF.LVM, CI = CI, evidence = evidence))
  }

}
