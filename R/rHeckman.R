#' Data generation from the Heckman Selection model (Normal, Student-t or CN)
#'
#' `rHeckman()` generates a random sample from the Heckman selection model (Normal, Student-t or CN).
#'
#' @param x A covariate matrix for the response y.
#' @param w A covariate matrix for the missing indicator cc.
#' @param beta Values for the beta vector.
#' @param gamma Values for the gamma vector.
#' @param sigma2 Value for the variance.
#' @param rho Value for the dependence between the response and missing value.
#' @param nu When using the t- distribution, the initial value for the degrees of freedom.
#' When using the CN distribution, the initial values for the proportion of bad observations and the degree of contamination.
#' @param family The family to be used (Normal, T, or CN).
#' @return Return an object with the response (y) and missing values (cc).
#'
#' @examples
#' \donttest{
#' n <- 100
#' rho <- .6
#' cens <- 0.25
#' nu <- 4
#' set.seed(20200527)
#' w <- cbind(1,runif(n,-1,1),rnorm(n))
#' x <- cbind(w[,1:2])
#'
#' family <- "T"
#' c <- qt(cens, df=nu)
#'
#' sigma2 <- 1
#' beta <- c(1,0.5)
#' gamma<- c(1,0.3,-.5)
#' gamma[1] <- -c*sqrt(sigma2)
#'
#' data <- rHeckman(x,w,beta,gamma,sigma2,rho,nu,family=family)
#' }
#' @export
rHeckman<-function(x,w,beta,gamma,sigma2,rho,nu=4,family="T"){
n <- nrow(x)
rhoa <- rho*sqrt(sigma2)

if(family=="Normal" || family=="normal"){
Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
#' @importFrom mvtnorm rmvnorm
errorTerms<- mvtnorm::rmvnorm(n, c( 0, 0 ), Sigma)
resp<- cbind(x%*%beta,w%*%gamma)+ errorTerms
}

if(family=="T" || family=="t"){
Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
#' @importFrom mvtnorm rmvt
errorTerms<- mvtnorm::rmvt(n, Sigma,df=nu)
resp<- cbind(x%*%beta,w%*%gamma)+ errorTerms
}

if(family=="CN" || family=="cn"){# Contaminated Normal
  cond1 = nu[1]<0 || nu[1]>1
  cond2 = nu[2]<0 || nu[2]>1
  if (length(nu) != 2) {
    stop("initial vector of length 2 must be provided for nu when using the CN distribution!")
  }else if (cond1 || cond2){
    stop("both components of the vector nu must be between 0 and 1 when using the CN distribution!")
  }
  #' @importFrom stats runif
  p <- runif(n)
  u <- rep(1,n)
  u[p<nu[1]] <- nu[2]
  Sigma<- matrix(c(sigma2, rhoa, rhoa, 1 ), ncol = 2)
  errorTerms<- rmvnorm(n, c( 0, 0 ), Sigma)/sqrt(u)
  resp<- cbind(x%*%beta,w%*%gamma)+ errorTerms
}



cc<-(resp[,2]> 0)+0
resp[cc==0,1]<-0

return=list(y=resp[,1],cc=cc)
}

