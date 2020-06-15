#' Title
#'
#' @param Data sdf
#' @param Model sdf
#' @param Prior sdf
#' @param R sdf
#' @param burnin sdf
#'
#' @return sdf
#' @export
#'
BFMediate = function(Data, Model, Prior, R, burnin){  # BF){
  # library(bayesm)
  # source('PartialMed_new.R', echo = F)
  # source('Mediation_Ordered_Multi_Merr.R', echo = F)
  # source('Measurement_Multi_Moment_Gen_New.R', echo = F)
  # source('rordprobitGibbs_me_M_multi_merr.R', echo = F)
  # source('rordprobitGibbs_me_multi_merr.R', echo = F)
  # source('Mediation_Ordered_Multi_Merr.R', echo = F)
  # source('MargeLik_SD.R')

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
  M = Data$M
  Y = Data$Y
  N = length(Data$X);
  evidence = NULL

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

  # B&K test results
  std_res = ifelse(BK2$coefficients["X","Pr(>|t|)"]<.05,"Reject", "Fail to reject")

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

  # BF simple
  if(Model == 'Simple'){
    out.simple = PartialMed(Data = Data, R=R,pars = list(A_M=A_M, A_Y=A_Y))
    beta_1 = out.simple$beta_1[,2]
    beta_2 = out.simple$beta_2[,2]
    beta_3 = out.simple$beta_2[,3]
    # BF.simple = exp(BFSD(Data = Data, post = out.simple , model = 'cont', prior = list(A=A_Y[3]),burnin = 0))
    BF.simple = exp(BFSD(post = out.simple , prior = A_Y[3], burnin = burnin))
    Bayes.CI.Indirect = round(as.vector(quantile(beta_1*beta_2,probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(beta_3,probs = c(.025,.975))),2)

    if(BF.simple>1) evidence = ifelse(BF.simple>100,"Decisive",
                                      ifelse(BF.simple>10,"strong",
                                             ifelse(BF.simple>3.2,"Substantial","Not worth more than a bare mention")))

    if(BF.simple<1) evidence = ifelse(1/BF.simple>100,"Decisive",
                                      ifelse(1/BF.simple>10,"strong",
                                             ifelse(1/BF.simple>3.2,"Substantial","Not worth more than a bare mention")))

    CI = as.character(ifelse((quantile(beta_3,probs = .025)>0) | (quantile(beta_3,probs = .975)<0),"Reject","Accept"))

    return(list(BK.beta_1 = BK.beta_1, BK.beta_1.sd = BK.beta_1.sd, BK.beta_2 = BK.beta_2, BK.beta_2.sd = BK.beta_2.sd, BK.beta_3 = BK.beta_3, BK.beta_3.sd = BK.beta_3.sd, BK.beta_3.pvalue = BK.beta_3.pvalue, std_res = std_res, evidence = evidence,
                PH.mean.Indirect = PH.mean.Indirect, PH.CI.Indirect = PH.CI.Indirect, PH.CI.Direct = PH.CI.Direct, Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
                beta_1 = beta_1, beta_2 = beta_2, beta_3 = beta_3, BF.simple = BF.simple, CI = CI))
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

    # stanfit = stan(file="Measurement_Multi.stan",
    #                data = list(n=N,M_ind = m_ind, Y_ind = y_ind,A_M=sqrt(A_M), A_Y=sqrt(A_Y), X = X, m_star=t(m_star),y_star=t(y_star)),
    #                pars = c("beta_1","beta_2","beta_3","ssq_M","ssq_Y","ssq_m_star","ssq_y_star","beta_0_Y","beta_0_M","lambda","tau","M","Y"),
    #                chains = 1,iter = R, warmup = burnin)
    #
    # out.multi = rstan::extract(stanfit)
    # moments = Measurement_Multi_Moment_Gen_New(Data = Data, post = out.multi, pars = list(A=1/A_Y), R = R - burnin)
    out.muti = MeasurementCont(Data = Data, Prior = list(A_M = A_M, A_Y = A_Y),R=R, burnin = burnin)
    BF.LVM = exp(BFSD(post = out.multi , prior = A_Y[3], burnin = 0)) #we already accounted for burnin in estimation
    Bayes.CI.Indirect = round(as.vector(quantile(out.multi$beta_1*out.multi$beta_2,probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(out.multi$beta_3,probs = c(.025,.975))),2)

    if(BF.LVM>1) evidence = ifelse(BF.LVM>100,"Decisive",
                                     ifelse(BF.LVM>10,"Strong",
                                            ifelse(BF.LVM>3.2,"Substantial","Not worth more than a bare mention")))

    if(BF.LVM<1) evidence = ifelse(1/BF.LVM>100,"Decisive",
                                     ifelse(1/BF.LVM>10,"Strong",
                                            ifelse(1/BF.LVM>3.2,"Substantial","Not worth more than a bare mention")))

    CI = as.character(ifelse((quantile(out.multi$beta_3, probs = .025)>0) | (quantile(out.multi$beta_3,probs = .975)<0),"Reject","Accept"))

    return(list(BK.beta_1 = BK.beta_1, BK.beta_1.sd = BK.beta_1.sd, BK.beta_2 = BK.beta_2, BK.beta_2.sd = BK.beta_2.sd, BK.beta_3 = BK.beta_3, BK.beta_3.sd = BK.beta_3.sd, BK.beta_3.pvalue = BK.beta_3.pvalue, std_res = std_res, evidence = evidence,
                PH.mean.Indirect = PH.mean.Indirect, PH.CI.Indirect = PH.CI.Indirect, PH.CI.Direct = PH.CI.Direct, Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
                out.multi = out.multi, BF.LVM = BF.LVM, CI = CI))
  }

  ################################
  #### Categorical Data Multi ####
  ################################

  # BF only M categorical
  if(Model == 'MCat'){

    Mcut = max(Data$m_star) +1
    Data_cat=list(X=cbind(rep(1,length(Data$X)),Data$X), m_star=as.matrix(Data$m_star), Y= as.matrix(Data$Y) ,k=Mcut-1, M_ind=dim(Data$m_star)[2])
    out.multi = MeasurementMCat(Data=Data_cat, Prior = Prior, R=R) #rordprobitGibbs_me_M_multi_merr_cpp(Data=Data_cat, Mcmc=Mcmc)
    BF.LVM = exp(BFSD(post = out.multi , prior = A_Y[3], burnin = burnin))
    Bayes.CI.Indirect = round(as.vector(quantile(out.multi$betadraw[,2]*out.multi$beta_2_draw[,2],probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(out.multi$beta_2_draw[,3],probs = c(.025,.975))),2)

    if(BF.LVM>1) evidence = ifelse(BF.LVM>100,"Decisive",
                                     ifelse(BF.LVM>10,"Strong",
                                            ifelse(BF.LVM>3.2,"Substantial","Not worth more than a bare mention")))

    if(BF.LVM<1) evidence = ifelse(1/BF.LVM>100,"Decisive",
                                     ifelse(1/BF.LVM>10,"Strong",
                                            ifelse(1/BF.LVM>3.2,"Substantial","Not worth more than a bare mention")))

    CI = as.character(ifelse((quantile(out.multi$beta_2_draw[,3], probs = .025)>0) | (quantile(out.multi$beta_2_draw[,3],probs = .975)<0),"Reject","Accept"))

    return(list(BK.beta_1 = BK.beta_1, BK.beta_1.sd = BK.beta_1.sd, BK.beta_2 = BK.beta_2, BK.beta_2.sd = BK.beta_2.sd, BK.beta_3 = BK.beta_3, BK.beta_3.sd = BK.beta_3.sd, BK.beta_3.pvalue = BK.beta_3.pvalue, std_res = std_res, evidence = evidence,
                PH.mean.Indirect = PH.mean.Indirect, PH.CI.Indirect = PH.CI.Indirect, PH.CI.Direct = PH.CI.Direct, Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
                out.multi = out.multi, BF.LVM = BF.LVM, CI = CI))
  }
  # BF only Y categorical
  if(Model == 'YCat'){


    Ycut = max(as.matrix(Data$y_star)[,1]) +1
    Data_cat=list(X=cbind(rep(1,length(Data$X)),Data$M,Data$X), y = as.matrix(Data$y_star) ,k=Ycut-1, Y_ind=dim(as.matrix(Data$y_star))[2])
    Mcmc=list(R=R)
    out.multi = MeasurementYCat(Data=DataMeasurementYCat, Prior=Prior, R=10000) #rordprobitGibbs_me_multi_merr_cpp(Data=Data_cat, Mcmc=Mcmc)
    BF.LVM = BF.LVM = exp(BFSD(post = out.multi , prior = A_Y[3], burnin = burnin))
    #computing beta_1 posterior using runiregGibbs to compute indirect effect Bayesian CI
    # outMX = runiregGibbs(Data=list(y=Data$M, X=cbind(rep(1,length(Data$X),Data$X))), Mcmc=list(R=R))
    Bayes.CI.Indirect = round(as.vector(quantile(out_MX$betadraw[,2]*out.multi$betadraw[,2], probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(out.multi$betadraw[,3],probs = c(.025,.975))),2)

    if(BF.LVM>1) evidence = ifelse(BF.LVM>100,"Decisive",
                                     ifelse(BF.LVM>10,"Strong",
                                            ifelse(BF.LVM>3.2,"Substantial","Not worth more than a bare mention")))

    if(BF.LVM<1) evidence = ifelse(1/BF.LVM>100,"Decisive",
                                     ifelse(1/BF.LVM>10,"Strong",
                                            ifelse(1/BF.LVM>3.2,"Substantial","Not worth more than a bare mention")))

    CI = as.character(ifelse((quantile(out.multi$betadraw[,3], probs = .025)>0) | (quantile(out.multi$betadraw[,3],probs = .975)<0),"Reject","Accept"))

    return(list(BK.beta_1 = BK.beta_1, BK.beta_1.sd = BK.beta_1.sd, BK.beta_2 = BK.beta_2, BK.beta_2.sd = BK.beta_2.sd, BK.beta_3 = BK.beta_3, BK.beta_3.sd = BK.beta_3.sd, BK.beta_3.pvalue = BK.beta_3.pvalue, std_res = std_res, evidence = evidence,
                PH.mean.Indirect = PH.mean.Indirect, PH.CI.Indirect = PH.CI.Indirect, PH.CI.Direct = PH.CI.Direct, Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
                out.multi = out.multi, BF.LVM = BF.LVM, CI = CI))
  }

  # BF both M and Y categorical
  if(Model == 'MYCat'){

    Mcut = max(Data$m_star) +1
    Ycut = max(Data$y_star) +1
    Data_cat=list(X=cbind(rep(1,length(Data$X)),Data$X), m_star=as.matrix(Data$m_star), y_star=as.matrix(Data$y_star), k_M = Mcut-1, k_Y=Ycut-1, M_ind=dim(as.matrix(Data$m_star))[2], Y_ind=dim(as.matrix(Data$y_star))[2])
    out.multi = MeasurementYCat(Data=DataMeasurementYCat, Prior=Prior, R=10000) #Mediation_Ordered_Multi_Merr(Data=Data_cat, Mcmc=Mcmc)
    BF.LVM = exp(BFSD(post = out.multi , prior = A_Y[3], burnin = burnin))
    Bayes.CI.Indirect = round(as.vector(quantile(out.multi$betadraw[,2]*out.multi$beta_2_draw[,2],probs = c(.025,.975))),2)
    Bayes.CI.Direct = round(as.vector(quantile(out.multi$beta_2_draw[,3],probs = c(.025,.975))),2)

    if(BF.LVM>1) evidence = ifelse(BF.LVM>100,"Decisive",
                                     ifelse(BF.LVM>10,"Strong",
                                            ifelse(BF.LVM>3.2,"Substantial","Not worth more than a bare mention")))

    if(BF.LVM<1) evidence = ifelse(1/BF.LVM>100,"Decisive",
                                     ifelse(1/BF.LVM>10,"Strong",
                                            ifelse(1/BF.LVM>3.2,"Substantial","Not worth more than a bare mention")))

    CI = as.character(ifelse((quantile(out.multi$beta_2_draw[,3], probs = .025)>0) | (quantile(out.multi$beta_2_draw[,3],probs = .975)<0),"Reject","Accept"))

    return(list(BK.beta_1 = BK.beta_1, BK.beta_1.sd = BK.beta_1.sd, BK.beta_2 = BK.beta_2, BK.beta_2.sd = BK.beta_2.sd, BK.beta_3 = BK.beta_3, BK.beta_3.sd = BK.beta_3.sd, BK.beta_3.pvalue = BK.beta_3.pvalue, std_res = std_res, evidence = evidence,
                PH.mean.Indirect = PH.mean.Indirect, PH.CI.Indirect = PH.CI.Indirect, PH.CI.Direct = PH.CI.Direct, Bayes.CI.Indirect = Bayes.CI.Indirect, Bayes.CI.Direct = Bayes.CI.Direct,
                out.multi = out.multi, BF.LVM = BF.LVM, CI = CI))
  }


}


#Partialmed_new -> PartialMed
#Margelik_SD --> BFSD,  model : measurement_sd --> cont, ordered--> cat
##change A_M and A_Y from prior SD to prior variance:
#partialmed_new
#measurement_multi.stan
