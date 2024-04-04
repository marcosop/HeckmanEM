#' Envelope for the Heckman Selection model
#'
#' `HeckmanEM.envelope()` plots the envelope for the fitted Heckman selection model.
#'
#' @param obj An object of the class HeckmanEM.
#' @param envelope The envelope coverage percentage.
#' @param ... Other option for chart.QQPlot from PerformanceAnalytics package.
#' @return A residual plot of the fitted data and its envelope.
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
#'                  im = TRUE, criteria = TRUE)
#' HeckmanEM.envelope(res, ylab="Normalized Quantile Residuals",xlab="Standard normal quantile",
#'                    line="quartile", col=c(20,1), pch=19, ylim = c(-5,4))
#' }
#' @export
HeckmanEM.envelope <- function(obj, envelope=0.95, ...)
{

  if(!inherits(obj,"HeckmanEM")) stop("Only \"HeckmanEM\" objects accepted!")
  if (obj$family != "Normal" && obj$family !="normal" && obj$family !="T" && obj$family !="t"  && obj$family !="CN" && obj$family !="cn") stop("Family not recognized")

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

  res <- resGCS(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, family)
  #res <- resMT(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, family) #old martingale envelopes
  if(family=="Normal" || family=="normal")
  {
    #' @importFrom PerformanceAnalytics chart.QQPlot
    PerformanceAnalytics::chart.QQPlot(res$resmt,distribution='norm',envelope=envelope,main="SLn", ...)
  }
  if(family=="T" || family=="t")
  {
    #' @importFrom PerformanceAnalytics chart.QQPlot
    PerformanceAnalytics::chart.QQPlot(res$resmt,distribution='norm',envelope=envelope,main="SLt", ...)
  }
  if(family=="CN" || family=="cn")
  {
    #' @importFrom PerformanceAnalytics chart.QQPlot
    PerformanceAnalytics::chart.QQPlot(res$resmt,distribution='norm',envelope=envelope,main="SLcn", ...)
  }
}

##-------------------------------##
## Martingale type residual (MT) ##
##-------------------------------##
resMT <- function(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, family)
{

  if(family == "normal" || family == "Normal") type <- "Normal"

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
    resmt[i] <- sqrt(-2*((1-cc[i])*log(1-cc[i]-resm[i])+resm[i]))*sign(resm[i])
  }

  return(list(resm=resm , resmt=resmt))
}

## ---------------------------------------------------- ##
## Envelopes of the normalized qualtile residuals       ##
## ---------------------------------------------------- ##
resGCS <- function(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, family= "Normal")
{

  if(family == "normal" || family == "Normal") type <- "Normal"
  if(family == "T" || family == "t") type <- "T"
  if(family == "cn" || family == "CN") type <- "CN"

  n<-nrow(x)
  y<-matrix(y,n,1)
  p<-ncol(x)
  q<-ncol(w)

  sigma2 <- sigma^2
  rhoa<- sigma*rho

  resm <- resmt <- c()

  S <-E<-matrix(0,n,1)
  if(type=="Normal")
  {
    mu1 <- x%*%beta #
    mu2<- w%*%gama
    mu12.1<- mu2+rho/sigma*(y-mu1)
    sigma12.1<- 1-rho^2
    Sigma<-matrix(c(sigma^2,sigma*rho,sigma*rho,1),2,2)

    for(i in 1:n){
      #' @importFrom MomTrunc pmvESN
      aux<-MomTrunc::pmvESN(c(-Inf,0),c(y[i],Inf),c(mu1[i],mu2[i]),Sigma,c(0,0),0)*cc[i]+pnorm(y[i]-mu2[i])*(1-cc[i])
      S[i] <- stats::qnorm(aux)
      E[i]  <-  -log(1-aux)
    }
  }
  if(type=="T")
  {
    mu1 <- x%*%beta # - sqrt(2/pi)*Delta
    mu2<- w%*%gama
    mu12.1<- mu2+rho/sigma*(y-mu1)
    sigma12.1<- 1-rho^2
    Sigma<-matrix(c(sigma^2,sigma*rho,sigma*rho,1),2,2)

    for(i in 1:n){
      #' @importFrom MomTrunc pmvEST
      aux<-MomTrunc::pmvEST(c(-Inf,0),c(y[i],Inf),c(mu1[i],mu2[i]),Sigma,c(0,0),0,nu)*cc[i]+pt(y[i]-mu2[i],nu)*(1-cc[i])
      S[i] <- stats::qnorm(aux)
      E[i]  <-  -log(1-aux)
    }
  }

  if(type=="CN"){

    nu1<-nu[1]
    nu2<-nu[2]
    mu1       <- x%*%beta # - sqrt(2/pi)*Delta
    mu2       <- w%*%gama
    mu12.1    <- mu2 + rho/sigma*(y - mu1)
    sigma12.1 <- 1 - rho^2
    Sigma     <- matrix(c(sigma^2, sigma*rho, sigma*rho, 1), 2, 2)

    for(i in 1:n){
      #' @importFrom MomTrunc pmvEST
      aux.1<- MomTrunc::pmvESN(c(-Inf, 0), c(y[i], Inf), c(mu1[i], mu2[i]), Sigma/nu2, c(0, 0), 0)
      aux.2<- MomTrunc::pmvESN(c(-Inf, 0), c(y[i], Inf), c(mu1[i], mu2[i]), Sigma, c(0, 0), 0)
      aux  <- (nu1*aux.1+(1-nu1)*aux.2)*cc[i] + (nu1* pnorm(y[i], mu2[i],sqrt(1/nu2))+(1-nu1)* pnorm(y[i] - mu2[i]))*(1 - cc[i])
      S[i] <- stats::qnorm(aux)
      E[i] <- -log(1 - aux)
    }
  }



  resmt<-S
  resm <-E

  return(list(resm=resm , resmt=resmt))
}


