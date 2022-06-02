#  
#'  Dataset of Sussman & O'Brien (2017) Study 3 
#'
#'  Contains the data collected to investigate whether the effect of saving intention (e.g. for children's education vs. for a vacation) on 
#'  the desire to withdraw from the saving account is mediated by a sense of personal responsibility.
#'
#' @format A dataframe with 408 observations and 6 variables:
#' 
#' \describe{
#' \item{...$CondGroup}{ Experimentally manipulated type of saving (X) }
#' \item{...$Irresponsible1}{ Perceived responsibility indicator 1 (m*_1) }
#' \item{...$Irresponsible2}{ Perceived responsibility indicator 2 (m*_2) }
#' \item{...$Irresponsible3}{ Perceived responsibility indicator 3 (m*_3) }
#' \item{...$IrresponsibleScale}{ Perceived responsibility composite measure }
#' \item{...$value_log}{ Compensation demanded (Y) }
#' }
#' 
#' @references {Sussman, A. B., & O'brien, R. L. (2016). Knowing when to spend: Unintended financial consequences of earmarking to encourage savings. Journal of Marketing Research, 53(5), 790-803.} 
#'
#' @examples 
#' 
#' data(Sussman_Obrien)
#' 
#' Data = NULL
#' Data$X = Sussman_Obrien$CondGroup
#' Data$M = Sussman_Obrien$IrresponsibleScale
#' Data$Y = scale(Sussman_Obrien$BorrowAmount)  
#' Data$m_tilde = cbind(Sussman_Obrien$Irresponsible1, Sussman_Obrien$Irresponsible2, Sussman_Obrien$Irresponsibility3_R)
#' 
#' # Setting priors
#' A_M = c(100,100)
#' A_Y = c(100,100,1)
#' 
#' # Computing Bayes factor for the simple mediation model using the composite measure of the mediator
#' out = Mediate(Data = Data, Model = 'Simple', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' 
#' # Parameter estimates
#' colMeans(out$Simple$beta_M)
#' colMeans(out$Simple$beta_Y)
#' out$Simple$Indirect_CI
#' out$Simple$Direct_CI
#' 
#' # Bayes factor
#' out$Simple$BF
#' out$Simple$evidence
#' 
#' # Computing Bayes factor for the latent variable mediation model using the indicators of the mediator
#' out_lvm = Mediate(Data = Data, Model = 'MCat', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' 
#' # Parameter estimates 
#' colMeans(out_lvm$beta_M)
#' colMeans(out_lvm$beta_Y)
#' out_lvm$Indirect_CI
#' out_lvm$Direct_CI
#'
#' # Bayes factor
#' out_lvm$BF 
#' out_lvm$evidence
#' 
#' # Computing Bayes factor for the latent variable mediation model with reverse causal ordering X->Y->M
#' 
#' # specifying Y as M, and Y indicators as M indicators
#' Data$y_tilde = Data$m_tilde
#' Data$m_tilde = NULL
#' Data$M = Data$Y
#' Data$Y = NULL
#' 
#' out_lvm_rev =  Mediate(Data = Data, Model = 'YCat', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' 
#' # Bayes factor
#' out_lvm_rev$BF
#' out_lvm_rev$evidence
"Sussman_Obrien"
