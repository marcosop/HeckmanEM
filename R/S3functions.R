#' @export
print.HeckmanEM <- function(x, ...){
  cat("\n*Family:", x$family,"\n")
  cat("*sample size:", length(x$y),"\n")
  cat("*beta:", x$beta,"\n")
  cat("*gamma:", x$gamma,"\n")
  cat("*sigma:", x$sigma,"\n")
  cat("*rho:", x$rho,"\n")
  if(x$family == "T" || x$family == "t" || x$family == "CN" || x$family == "cn") cat("*nu:",x$nu,"\n")
  #print(x)
}


#' @export
summary.HeckmanEM <- function(object, ...){
  #object = aa
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
  cat('Parameter Estimates\n')
  cat('------------------\n')
  cat('\n')
  # key <- FALSE
  if(object$family == "Normal" || object$family == "normal") est=c(object$beta,object$gamma,object$sigma,object$rho)
  if (object$family == "T" || object$family == "t"){
    est=c(object$beta,object$gamma,object$sigma,object$rho,object$nu)
    #key = TRUE
  }
  if (object$family == "CN" || object$family == "cn"){
    est=c(object$beta,object$gamma,object$sigma,object$rho,object$nu)
    if (length(est)>length(object$sd)) key = TRUE
  }
  l = length(est)
  lab = c(paste('beta',0:(length(object$beta)-1),sep=''),paste('gamma',0:(length(object$gamma)-1),sep=''),"sigma","rho")
  if(object$family == "T" || object$family == "t") lab <- c(lab,"nu")
  if(object$family == "CN" || object$family == "cn") lab <- c(lab,"nu1", "nu2")
  tab <- matrix(NA, ncol=l, nrow=2)
  tab[1,] = t(round(est,digits))
  if (object$family == "T" || length(est)>length(object$sd)){
    tab[2,] = round(c(object$sd,NA),digits)
  }else{
    tab[2,] = round(object$sd,digits)
  }
  colnames(tab)=t(lab)
  rownames(tab)=c("Estimate","Standard Error")
  print(tab)
}


