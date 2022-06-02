#  
#'  Dataset of Tonietto & Malkoc (2016) Study 4
#'
#'  Contains the data collected to investigate the serial mediation:
#'  Type of scheduling -> Free-flowing -> Work construal -> Anticipation utility (X->M1->M2->Y).
#'
#' @format A dataframe with 141 observations and 9 variables:
#' 
#' \describe{
#' \item{...$Condition}{ Experimentally manipulated type of scheduling (X) }
#' \item{...$Freeflowing}{ Free flowing indicator 1 (m_1*_1) }
#' \item{...$Flexible}{ Free flowing indicator 2 (m_1*_2) }
#' \item{...$Scale_Freeflowing}{ Free flowing composite measure }
#' \item{...$Chore}{ Work construal indicator 1 (m_2*_1) }
#' \item{...$Effortful}{ Work construal indicator 2 (m_2*_2) }
#' \item{...$WorkLike}{ Work construal indicator 3 (m_2*_3) }
#' \item{...$Scale_Work}{ Work construal composite measure }
#' \item{...$Scale_Anticipation_Utility}{  Anticipation utility composite measure }
#' }
#' 
#' @references {Tonietto, G. N., & Malkoc, S. A. (2016). The calendar mindset: Scheduling takes the fun out and puts the work in. Journal of Marketing Research, 53(6), 922-936.} 
#'
#' @examples 
#' 
#' data(Tonietto_Malkoc)
#' 
#' # Setting priors
#' A_M = c(100,100)
#' A_Y = c(100,100,1)
#' 
#' # Computing Bayes factor for each of the simple mediation chains:
#' # X->M1->M2, X->M1->Y, X->M2->Y, and M1->M2->Y 
#' # Substantial evidence in favor of conditional independence for all the above chains is evidence in favor of 
#' # the serial mediation X->M1->M2->Y 
#' 
#' # X->M1->M2
#' 
#' Data = NULL
#' Data$X = ifelse(Tonietto_Malkoc$Condition=="Impromptu",2,-1)    # C1 contrast coding of X (based on the paper)
#' Data$M = Tonietto_Malkoc$Scale_Freeflowing
#' Data$Y = Tonietto_Malkoc$Scale_Work
#' 
#' out1 = Mediate(Data = Data, Model = 'Simple', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' out1$Simple$BF
#' out1$Simple$evidence
#' 
#' # X->M1->Y
#' 
#' Data = NULL
#' Data$X = ifelse(Tonietto_Malkoc$Condition=="Impromptu",2,-1)
#' Data$M = Tonietto_Malkoc$Scale_Freeflowing
#' Data$Y = Tonietto_Malkoc$Scale_Anticipation_Utility
#' 
#' out2 = Mediate(Data = Data, Model = 'Simple', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' out2$Simple$BF
#' out2$Simple$evidence
#' 
#' # X->M2->Y
#' 
#' Data = NULL
#' Data$X = ifelse(Tonietto_Malkoc$Condition=="Impromptu",2,-1)
#' Data$M = Tonietto_Malkoc$Scale_Work
#' Data$Y = Tonietto_Malkoc$Scale_Anticipation_Utility
#' 
#' out3 = Mediate(Data = Data, Model = 'Simple', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' out3$Simple$BF
#' out3$Simple$evidence
#' 
#' # M1->M2->Y
#' 
#' Data = NULL
#' Data$X = Tonietto_Malkoc$Scale_Freeflowing
#' Data$M = Tonietto_Malkoc$Scale_Work
#' Data$Y = Tonietto_Malkoc$Scale_Anticipation_Utility
#' 
#' out4 = Mediate(Data = Data, Model = 'Simple', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' out4$Simple$BF
#' out4$Simple$evidence
#' 
#' # BF shows evidence against conditional independence in the chain M1->M2->Y (the results are robust to a latent variable model specification using the indicators)
#' # therefore we cannot provide evidence in favor of serial mediation. 
#' 
#' # Next we explore a different data-generating process where M1 and M2 are indicators of a latent mediator
#'
#' Data$X = ifelse(Tonietto_Malkoc$Condition=="Impromptu",2,-1)
#' # using the indicators of (reverse) Freeflowing and work construal as indicators for a common latent mediator
#' Data$m_tilde = cbind(Tonietto_Malkoc$Chore,Tonietto_Malkoc$Effortful,Tonietto_Malkoc$WorkLike,8-Tonietto_Malkoc$Flexible,8-Tonietto_Malkoc$Freeflowing)
#' Data$Y = Tonietto_Malkoc$Scale_Anticipation_Utility
#' 
#' # Computing Bayes factor for the new latent variable model
#' outnew = Mediate(Data = Data, Model = 'MCat', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' outnew$BF
#' outnew$evidence
"Tonietto_Malkoc"