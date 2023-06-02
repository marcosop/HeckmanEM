#' Model selection criteria for the Heckman Selection model
#'
#' `HeckmanEM.criteria()` calculates the AIC, AICc, BIC selection criteria for the fitted Heckman selection model.
#'
#' @param obj An object of the class HeckmanEM.
#' @return The calculated AIC, AICc, and BIC for the parameters of the fitted model.
#'
#' @examples
#' \donttest{
#' n <- 100
#' family <- "T"
#' nu <- 4
#' rho <- .6
#' cens <- .25
#'
#' set.seed(20200527)
#' w <- cbind(1,runif(n,-1,1),rnorm(n))
#' x <- cbind(w[,1:2])
#' c <- qt(cens, df=nu)
#'
#' sigma2 <- 1
#'
#' beta <- c(1,0.5)
#' gamma <- c(1,0.3,-.5)
#' gamma[1] <- -c*sqrt(sigma2)
#'
#' set.seed(1)
#' datas <- rHeckman(x,w,beta,gamma,sigma2,rho,nu,family=family)
#' y <- datas$y
#' cc <- datas$cc
#'
#' res <- HeckmanEM(y, x, w, cc, nu = 4, family = "Normal", error = 1e-05, iter.max = 500,
#'                  im = TRUE, criteria = FALSE)
#' cr <- HeckmanEM.criteria(res)
#' }
#' @export
HeckmanEM.criteria <- function(obj){

  if(!inherits(obj,"HeckmanEM")) stop("Only \"HeckmanEM\" objects accepted!")
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
