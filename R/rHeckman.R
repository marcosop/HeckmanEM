#' Generate data from the Heckman Selectio model (Normal or Student-t).
#'
#' `rHeckman()` generate a ramdom sample from the Heckman selection model (Normal or Student-t).
#'
#' @param x A covariate matrix for the response y.
#' @param w A covariate matrix for the missing indicator cc.
#' @param beta Values for the beta vector.
#' @param gamma Values for the gamma vector.
#' @param sigma2 Value for the variance.
#' @param rho Value for the dependence between the reponse and missing value.
#' @param nu Value for the degrees of freedom.
#' @param family The family to be used (Normal or T).
#' @return Return an object with the response (y) and missing values (cc).
#'
#' @examples
#' \dontrun{
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
#' set.seed(iter)
#' data <- rHeckman(x,w,beta,gamma,sigma2,rho,nu,family=family)
#' }
#' @export
rHeckman<-function(x,w,beta,gamma,sigma2,rho,nu=4,family="T"){
n<-nrow(x)
rhoa<- rho*sqrt(sigma2)

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

cc<-(resp[,2]> 0)+0
resp[cc==0,1]<-0

return=list(y=resp[,1],cc=cc)
}

