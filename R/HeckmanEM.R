#' Fit the Normal, Student-t or Contaminated normal Heckman Selection model
#'
#' `HeckmanEM()` fits the Heckman selection model.
#'
#' @param y A response vector.
#' @param x A covariate matrix for the response y.
#' @param w A covariate matrix for the missing indicator cc.
#' @param cc A missing indicator vector (1=observed, 0=missing) .
#' @param nu When using the t- distribution, the initial value for the degrees of freedom.
#' When using the CN distribution, the initial values for the proportion of bad observations and the degree of contamination.
#' @param family The family to be used (Normal, T or CN).
#' @param error The absolute convergence error for the EM stopping rule.
#' @param iter.max The maximum number of iterations for the EM algorithm.
#' @param im TRUE/FALSE, boolean to decide if the standard errors of the parameters should be computed.
#' @param criteria TRUE/FALSE, boolean to decide if the model selection criteria should be computed.
#' @param verbose TRUE/FALSE, boolean to decide if the progress should be printed in the screen.
#' @return An object of the class HeckmanEM with all the outputs provided from the function.
#'
#' @examples
#' n    <- 100
#' nu   <- 3
#' cens <- 0.25
#'
#' set.seed(13)
#' w <- cbind(1,runif(n,-1,1),rnorm(n))
#' x <- cbind(w[,1:2])
#' c <- qt(cens, df=nu)
#'
#' sigma2   <- 1
#' beta     <- c(1,0.5)
#' gamma    <- c(1,0.3,-.5)
#' gamma[1] <- -c*sqrt(sigma2)
#'
#' set.seed(1)
#' datas <- rHeckman(x,w,beta,gamma,sigma2,rho = 0.6,nu,family="T")
#' y <- datas$y
#' cc <- datas$cc
#' \donttest{
#' # Normal EM
#' res.N <- HeckmanEM(y, x, w, cc, family="Normal",iter.max = 50)
#' # Student-t EM
#' res.T <- HeckmanEM(y, x, w, cc, nu = 4, family="T", iter.max = 50)
#' }
#' @export
HeckmanEM <- function(y, x, w, cc, nu = 4, family="T", error = 1e-05,iter.max = 500, im=TRUE, criteria = TRUE, verbose = TRUE){

 if (family != "Normal" && family !="normal" && family !="T" && family !="t" && family !="CN" && family !="cn" ) stop("Family not recognized! Obly families allowed are: \"Normal\", \"T\" and \"CN\".")
 if(!is.vector(y)) stop("y must be a vector!")
 if(!is.vector(cc)) stop("y must be a vector!")

 if(is.vector(x)) x <- as.matrix(x)
 if(is.vector(w)) w <- as.matrix(w)
 if(!is.matrix(x)) stop("y must be a matrix!")
 if(!is.matrix(w)) stop("y must be a matrix!")

 if(family == "Normal" || family == "normal"){
   out <- EMn.alg(y=y, x=x, w=w, cc=cc, error = error, iter.max = iter.max, im=im, criteria = criteria, verbose=verbose)
 }

 if((family == "T" || family == "t")){
   if (length(nu) == 0) {
     stop("initial for nu must be provided!")
   } else{
     out <- EMt.alg(y=y, x=x, w=w, cc=cc, nu=nu, error = error, iter.max = iter.max, im=im, criteria = criteria, verbose=verbose)
     }
 }


 if((family == "CN" || family == "cn")){
   cond1 = nu[1]<0 || nu[1]>1
   cond2 = nu[2]<0 || nu[2]>1
   if (length(nu) != 2) {
     stop("initial vector of length 2 must be provided for nu when using the CN distribution!")
   }else if (cond1 || cond2){
     stop("both components of the vector nu must be between 0 and 1 when using the CN distribution!")
   }
   else{
     out <- EMcn.alg(y=y, x=x, w=w, cc=cc, nu=nu, error = error, iter.max = iter.max, im=im, criteria = criteria, verbose=verbose)
   }
 }

 return(out)
}
