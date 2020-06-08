#' Envelope for the Heckman Selectio model.
#'
#' `HeckmanEM.envelope()` plot the envelope for the fitted Heckman selection model.
#'
#' @param obj An object of the class HeckmanEM.
#' @param envelope The envelope corverage percertange.
#' @param ... Other option for chart.QQPlot from PerformanceAnalytics package.
#' @return A residual plot of the fitted data and its envelope.
#'
#' @examples
#' \dontrun{
#' res <- HeckmanEM(y, x, w, cc, nu = 4, family = "Normal", error = 1e-05, iter.max = 500, 
#'                  im = TRUE, criteria = TRUE)
#' HeckmanEM.envelope(res, ylab="Martingale-type residuals",xlab="Standard normal quantile",
#'                    line="quartile", col=c(20,1),pch=19)
#' }
#' @export
HeckmanEM.envelope <- function(obj, envelope=0.95, ...)
{

  if(class(obj) != "HeckmanEM") stop("Only \"HeckmanEM\" objects accepted!")
  if (obj$family != "Normal" && obj$family !="normal" && obj$family !="T" && obj$family !="t" ) stop("Family not recognized")

  y <- obj$y
  x <- obj$x
  w <- obj$w
  cc <- obj$cc
  beta <- obj$beta
  gama <- obj$gamma
  rho <- obj$rho
  sigma <- obj$sigma
  nu <- obj$nu
  family  <- obj$family

  res <- resMT(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, family)
  if(family=="Normal" || family=="normal")
  {
#' @importFrom PerformanceAnalytics chart.QQPlot
    PerformanceAnalytics::chart.QQPlot(res$resmt,distribution='norm',envelope=envelope,main="SLn", ...)
  }
  else
  {
    PerformanceAnalytics::chart.QQPlot(res$resmt,distribution='norm',envelope=envelope,main="SLt", ...)
  }
  
}

##-------------------------------##
## Martingale type residual (MT) ##
##-------------------------------##
resMT <- function(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, family)
{

  if(family == "normal" || family == "Normal") type <- "Normal"
  else type <- "T"
  n<-nrow(x)
  y<-matrix(y,n,1)
  p<-ncol(x)
  q<-ncol(w)
  cc<-1-cc

  sigma2 <- sigma^2
  rhoa<- sigma*rho

  resm <- resmt <- c()
  
  S <-matrix(0,n,1)
  if(type=="Normal")
  {
    mu1 <- x%*%beta 
    mu2<- w%*%gama
    mu12.1<- mu2+rho/sigma*(y-mu1)
    sigma12.1<- 1-rho^2
    S <- pnorm(-mu2)
  }
  if(type=="T")
  {
    mu1 <- x%*%beta 
    mu2<- w%*%gama
    mu12.1<- mu2+rho/sigma*(y-mu1)
    sigma12.1<- 1-rho^2
    S <- pt(-mu2,nu)
  }
  
  for (i in 1:n){
    resm[i] <-  1-cc[i]+log(S[i])
    resmt[i] <-   sqrt(-2*((1-cc[i])*log(1-cc[i]-resm[i])+resm[i]))*sign(resm[i])
  }
  
  return(list(resm=resm , resmt=resmt)) 
}

