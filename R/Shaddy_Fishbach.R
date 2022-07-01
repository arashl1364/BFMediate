#  
#'  Dataset of Shaddy & Fishbach (2017) Study 5 (Scenario: lose)
#'
#'  Contains the data collected to investigate whether the effect of losing a product in a bundle (vs. in isolation) on 
#'  the compensation demanded for the product loss is mediated by the need to replace that product
#'
#' @format A dataframe with 103 observations and 3 variables:
#' 
#' \describe{
#' \item{...$bundled_cond}{ Experimentally manipulated Bundling condition (X)}
#' \item{...$replacement}{ Perceived importance of replacement (M)}
#' \item{...$value_log}{ Compensation demanded (Y)}
#' }
#' 
#' @references {Shaddy, F., & Fishbach, A. (2017). Seller beware: how bundling affects valuation. Journal of Marketing Research, 54(5), 737-751.} 
#'
#' @examples 
#' 
#' data(Shaddy_Fishbach)
#' 
#' Data = NULL
#' Data$X = Shaddy_Fishbach$bundled_cond
#' Data$M = Shaddy_Fishbach$replacement
#' Data$Y = Shaddy_Fishbach$value_log
#' 
#' # Saving the dataset to use in the Shiny app (https://bfmediate.shinyapps.io/bfmediate_app/)
#' save(Data,file = "~/Shaddy_Fishbach.rda")    # the file path can be changed by replacing ~
#' 
#' # Setting priors
#' A_M = c(100,100)
#' A_Y = c(100,100,1)
#' 
#' # Computing Bayes factor for the simple mediation model
#' out = Mediate(Data = Data, Model = 'Simple', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' 
#' # Parameter estimates
#' colMeans(out$Simple$beta_M)
#' colMeans(out$Simple$beta_Y)
#' 
#' # Bayes factor
#' out$Simple$BF
#' out$Simple$evidence
#' 
#' # Computing the Bayes factor for the model with reverse MY causal direction (X->Y->M)
#' 
#' # Specifying Y as M and vice versa
#' temp = Data$M
#' Data$M = Data$Y
#' Data$Y = temp
#' 
#' out_rev = Mediate(Data = Data, Model = 'Simple', Prior = list(A_M = A_M, A_Y = A_Y),R=10000, burnin = 2000)
#' out_rev$Simple$BF
#' out_rev$Simple$evidence
"Shaddy_Fishbach"
