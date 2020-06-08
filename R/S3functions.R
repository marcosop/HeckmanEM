#' @export
print.HeckmanEM <- function(x, ...){
cat("\n*Family:", x$family,"\n")
cat("*sample size:", length(x$y),"\n")
cat("*beta:", x$beta,"\n")
cat("*gamma:", x$gamma,"\n")
cat("*sigma:", x$sigma,"\n")
cat("*rho:", x$rho,"\n")
if(x$family == "T" || x$family == "t") cat("*nu:", x$nu,"\n")
}

#' @export
summary.HeckmanEM <- function(object, ...){
var <- as.list(match.call())
if(!("digits" %in% names(var))) digits <- 4
else digits <- var$digits
cat('\n')
cat('---------------------------------------------------\n')
cat('  Heckman Selection Model \n')
cat('---------------------------------------------------\n')
cat('\n')
cat("*Family:", object$family,'\n')
cat("*sample size:", length(object$y),"\n")
cat('------------------\n')
cat('Paramter Estimates\n')
cat('------------------\n')
cat('\n')
key <- FALSE
if(object$family == "Normal" || object$family == "normal") est=c(object$beta,object$gamma,object$sigma,object$rho)
else{
 est=c(object$beta,object$gamma,object$sigma,object$rho,object$nu)
 key = TRUE
}
l = length(est)
lab = c(paste('beta',0:(length(object$beta)-1),sep=''),paste('gamma',0:(length(object$gamma)-1),sep=''),"sigma","rho")
if(object$family == "T" || object$family == "t") lab <- c(lab,"nu")
tab <- matrix(NA, ncol=l, nrow=2)
tab[1,] = t(round(est,digits))
tab[2,] = ifelse(key, round(c(object$sd,NA),digits), round(object$sd,digits)) 
colnames(tab)=t(lab)
rownames(tab)=c("Estimate","Standard Error")
print(tab)
}
