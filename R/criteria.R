#' Model selection criteria for the Heckman Selectio model
#'
#' `HeckmanEM.criteria()` calculate the AIC, AICc, BIC selection criteria for the fitted Heckman selection model.
#'
#' @param obj An object of the class HeckmanEM.
#' @return The calculated AIC, AICc, and BIC for the parameters of the fitted model.
#'
#' @examples
#' \dontrun{
#' res <- HeckmanEM(y, x, w, cc, nu = 4, family = "Normal", error = 1e-05, iter.max = 500, 
#'                  im = TRUE, criteria = FALSE)
#' cr <- HeckmanEM.criteria(res)
#' }
#' @export
HeckmanEM.criteria <- function(obj){

  if(class(obj) != "HeckmanEM") stop("Only \"HeckmanEM\" objects accepted!")
  if (obj$family != "Normal" && obj$family !="normal" && obj$family !="T" && obj$family !="t" ) stop("Family not recognized")

  n <- length(obj$y)
  p <- length(obj$beta)
  q <- length(obj$gamma)

  if(obj$family == "Normal" || obj$family == "normal") npar<-length(p+q+2)
  else npar<-length(p+q+3)

  lkante <- obj$logL
  ##Model comparison criteria
  AICc<- -2*lkante +2*npar
  AICcorr<- AICc + ((2*npar*(npar+1))/(n-npar-1))
  BICc <- -2*lkante +log(n)*npar

  return(list(AIC=AICc, AICc = AICcorr, BIC = BICc))
}
