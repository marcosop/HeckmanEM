#' Fit the Normal or Student-t Heckman Selectio model.
#'
#' `HeckmanEM()` fit the Heckman selection model.
#'
#' @param y A response vector.
#' @param x A covariate matrix for the response y.
#' @param w A covariate matrix for the missing indicator cc.
#' @param cc A missing incidator vector (1=obserced, 0=missing) .
#' @param nu The initial value for the degrees of freedom.
#' @param family The family to be used (Normal or T).
#' @param error The abslute convergence error for the EM stopping rule.
#' @param iter.max The maximum number of iterations for the EM algorithm.
#' @param im TRUE/FALSE, boolean to decide if the standard erros of the parameters should be computed.
#' @param criteria TRUE/FALSE, boolean to decide if the model selection criteria should be computed.
#' @return An object of the class HeckmanEM with all the outputs provided from the function.
#'
#' @examples
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
#' \donttest{
#' # Normal EM
#' res.N <- HeckmanEM(y, x, w, cc, nu = 4, family="Normal", error = 1e-05, iter.max = 50,
#'                    im=TRUE, criteria = TRUE)
#' # Student-t: EM
#' res.T <- HeckmanEM(y, x, w, cc, nu = 4, family="T", error = 1e-05, iter.max = 50,
#'                    im=TRUE, criteria = TRUE)
#' }
#' @export
HeckmanEM <- function(y, x, w, cc, nu = 4, family="T", error = 1e-05,iter.max = 500, im=TRUE, criteria = TRUE){

 if (family != "Normal" && family !="normal" && family !="T" && family !="t" ) stop("Family not recognized! Obly families allowed are: \"Normal\" and \"T\".")
 if(!is.vector(y)) stop("y must be a vector!")
 if(!is.vector(cc)) stop("y must be a vector!")

 if(is.vector(x)) x <- as.matrix(x)
 if(is.vector(w)) w <- as.matrix(w)
 if(!is.matrix(x)) stop("y must be a matrix!")
 if(!is.matrix(w)) stop("y must be a matrix!")

 if((family == "T" || family == "t") && length(nu) == 0) stop("initial for nu must be provided!")

 if(family == "Normal" || family == "normal"){
   out <- EMn.alg(y=y, x=x, w=w, cc=cc, error = error, iter.max = iter.max, im=im, criteria = criteria)
 }
 else out <- EMt.alg(y=y, x=x, w=w, cc=cc, nu=nu, error = error, iter.max = iter.max, im=im, criteria = criteria)

 return(out)
}
