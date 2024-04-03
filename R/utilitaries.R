################################################################################
## NOT SHARED FUNCTIONS
################################################################################
#' @importFrom stats as.formula dnorm dt optimize pnorm pt
################################################################################
### likelihood function  : Normal
################################################################################
likeL<-function(y, x, w, cc, beta, gama, Sigma){
sigma2<- Sigma[1,1]
sigma<- sqrt(sigma2)
rho<- Sigma[1,2]/sigma
med<-(y-x%*%beta)/sigma

fd <- ifelse(cc==1, stats::dnorm(med,log=T)+stats::pnorm((rho*med+w%*%gama)/sqrt(1-rho^2), log.p=T)-log(sigma), stats::pnorm(-w%*%gama,log.p=T))
return(sum(fd))
}

################################################################################
### likelihood function  : Student--t
################################################################################
likeLt<-function(nu, y, x, w, cc, beta, gama, Sigma){
sigma2<- Sigma[1,1]
sigma<- sqrt(sigma2)
rho<- Sigma[1,2]/sigma
med<-(y-x%*%beta)/sigma
aux1<-sqrt((nu+1)/(nu+med^2))*((rho*med+w%*%gama)/(sqrt(1-rho^2)))

fd<-ifelse(cc==1, stats::dt(med,nu,log=T)+stats::pt(aux1,nu+1,log.p=T)-log(sigma), stats::pt(-w%*%gama,nu,log.p=T))
return(sum(fd))
}

################################################################################
# generate initial values
################################################################################
init.values <- function(y,x,w,cc){

  ##generate initial values
  colnames(w) <- paste0("w",0:(ncol(w)-1))
  colnames(x) <- paste0("x",0:(ncol(x)-1))

  dat <- data.frame(y=y,cc=as.logical(cc),x,w)
  fy <- "y ~ -1 +"
  fy <- as.formula(paste(fy, paste(colnames(x),collapse="+"), sep=" "))

  fm <- "cc ~ -1 +"
  fm <- as.formula(paste(fm, paste(colnames(w),collapse="+"), sep=" "))

  # Two-step estimation
#' @importFrom sampleSelection heckit
  out <-  sampleSelection::heckit(fm, fy, data=dat, method = "2step" )

  beta <- out$coefficients[ncol(w)+1:ncol(x)]
  gamma <- out$coefficients[1:ncol(w)]
  sigma <- out$coefficients["sigma"]
  rho <- out$coefficients["rho"]
  return(list(beta=beta,gamma=gamma,sigma=sigma,rho=rho))
}

################################################################################
# Algoritmo Normal case: Using truncated moments
################################################################################
EMn.alg<-function(y, x, w, cc, error,iter.max, im, criteria,verbose){

n<-nrow(x)
y<-matrix(y,n,1)
p<-ncol(x)
q<-ncol(w)

ini <- init.values(y,x,w,cc)

beta <- ini$beta
gama <- ini$gamma
rho <- ifelse(abs(ini$rho) >= 1-1e-07, 0.5, ini$rho)
sigma <-  ifelse(ini$sigma <= 0, 1, ini$sigma)

sigma2<- sigma^2
rhoa<- sigma*rho
Sigma<- matrix(c(sigma2, rhoa, rhoa,1), ncol = 2)

betaC<-c(beta,gama)

 criterio <- 1
 count    <- 0
 lkante   <- likeL(y, x, w, cc, beta, gama, Sigma)


while((criterio > error) && (count <= iter.max)){
   count<-count+1
   if(verbose) cat("Iteration: ", count,"\r")

   mu1<- x%*%beta
   mu2<- w%*%gama


   suma1<- suma11<-suma12<-0 #matrix(0,2,2)
   suma2<- soma4<-matrix(0,p+q,p+q)
   suma3<- matrix(0,p+q,1)
 	for (i in 1:n){

       uy <-matrix(y[i],2,1)
       uyy<-matrix(0,2,2)

  		if(cc[i]==1){
				mu12.1<- mu2[i]+rho/sigma*(y[i]-mu1[i])
        sigma12.1<- 1-rho^2

         if(sigma12.1 == 0) sigma12.1 <- 0.0001
        #' @importFrom MomTrunc meanvarTMD
        MomNT<- MomTrunc::meanvarTMD(lower = 0,upper = Inf, mu = mu12.1, Sigma = sigma12.1, dist="normal")

        uy[2]<-MomNT$mean
        uyy[2,2]<-MomNT$varcov
        uyy<- uyy+uy%*%t(uy)

   		}
      else{

        MomNT1<- MomTrunc::meanvarTMD(lower = c(-Inf,-Inf),upper = c(Inf,0),mu = c(mu1[i],mu2[i]),Sigma = Sigma,dist="normal")
        uy <-MomNT1$mean
        uyy<-MomNT1$EYY

      }
     SigI <- solve(Sigma, tol=1e-30)
     aux41<- rbind(c(x[i,],0*w[i,]),c(x[i,]*0,w[i,]))
     aux5 <- aux41%*%betaC
     aux6<- uyy-aux5%*%t(uy)-t(aux5%*%t(uy))+aux5%*%t(aux5)
     suma1<- suma1+ aux6[1,1]-rhoa*(aux6[1,2]+aux6[2,1])+aux6[2,2]*rhoa^2
     suma11<-suma11+(aux6[1,2]+aux6[2,1])/2
     suma12<-suma12+aux6[2,2]
     suma2<- suma2+ t(aux41)%*%SigI%*%(aux41)
     suma3<- suma3+ t(aux41)%*%SigI%*%uy
  }

  suma2aux<-(suma2+t(suma2))/2

  betaC<- solve(suma2aux, tol=1e-30)%*%suma3
  psi<- suma1/n
  rhoa<- suma11/suma12
  sigma2<- psi+rhoa^2
  sigma<- sqrt(sigma2)
  rho<- rhoa/sigma
  Sigma<- matrix(c(sigma2,rhoa,rhoa,1),2,2)

  beta<- betaC[1:p]
  gama<- betaC[(p+1):(p+q)]


 lkante1<-lkante
 lkante<- likeL(y, x, w, cc, beta, gama, Sigma)
 criterio <- sqrt(abs(1-lkante1/lkante))

 if((rho >= 1 - 1e-07)&&(rho <= -1 + 1e-07)) {
     stop("EM did not coverge")#criterio<-0.0000000000001
  }
 }

 if(verbose) cat("\n")
 AICc<- NULL
 AICcorr <- NULL
 BICc <- NULL
 desvios <- NULL

 out <- list(y=y, x=x, w=w, cc=cc, beta=beta, gamma=gama, rho=rho, sigma=sigma, nu=NULL, sd=desvios, logL=lkante,AIC=AICc, AICc=AICcorr, BIC=BICc, family="Normal")
 class(out) <- "HeckmanEM"

 if(criteria){
  crit <- HeckmanEM.criteria(out)
  out$AIC <- crit$AIC
  out$AICc <- crit$AICc
  out$BIC <- crit$BIC
 }

 if(im){
   desvios <- HeckmanEM.infomat(out)
   out$sd <- desvios
 }

 return(out)
}


################################################################################
# Algorithm Student-t: Using truncated moments
################################################################################
EMt.alg<-function(y, x, w, cc, nu, error,iter.max,im,criteria,verbose){

#' @importFrom mvtnorm GenzBretz
GB = mvtnorm::GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)

n<-nrow(x)
y<-matrix(y,n,1)
p<-ncol(x)
q<-ncol(w)

ini <- init.values(y,x,w,cc)

beta <- ini$beta
gama <- ini$gamma
rho <- ifelse(abs(ini$rho) >= 1-1e-07, 0.5, ini$rho)
sigma <-  ifelse(ini$sigma <= 0, 1, ini$sigma)

sigma2<- sigma^2
rhoa<- sigma*rho
Sigma<- matrix(c(sigma2, rhoa, rhoa,1), ncol = 2)

betaC<-c(beta,gama)

 criterio <- 1
 count    <- 0
lkante   <- likeLt(nu, y, x, w, cc, beta, gama, Sigma)

while((criterio > error) && (count <= iter.max)){
   count<-count+1
   if(verbose) cat("Iteration: ", count,"\r")

   mu1<- x%*%beta
   mu2<- w%*%gama


   suma1<- suma11<-suma12<-0 #matrix(0,2,2)
   suma2<-soma4<-matrix(0,p+q,p+q)
   suma3<-matrix(0,p+q,1)

 	for (i in 1:n){

       uy <-matrix(y[i],2,1)
       uyy<-matrix(0,2,2)

  		if(cc[i]==1){
         PsiA<-Sigma*nu/(nu+2)
         nu1<-(nu+1)##X
         mu12.1<- mu2[i]+rho/sigma*(y[i]-mu1[i])
         sigma12.1<- 1-rho^2
         ScA<-nu/(nu+2)*sigma12.1

         Qy1<- (y[i]-mu1[i])^2/sigma2  #X
         Qy2<- (y[i]-mu1[i])^2/PsiA[1,1]  #X
          auxcte<-as.numeric((nu+Qy1)/(nu+1))
          auxcte1<-as.numeric((nu+2+Qy2)/(nu+3))

          Sc22<-auxcte*sigma12.1

          muUi<-mu12.1
          SigmaUi<-Sc22

          SigmaUiA<- auxcte1*ScA
          auxupper<- -muUi

          auxU1<- 1-pt(auxupper/sqrt(SigmaUiA),nu1+2)

          auxU2<- 1-pt(auxupper/sqrt(SigmaUi),nu1)



          MoMT<- MomTrunc::meanvarTMD(lower = 0,upper = Inf,mu = muUi, Sigma = SigmaUiA, dist="t",nu=nu1+2)


         vary <- matrix(0,2,2)
         vary[2,2]<- MomTrunc::meanvarTMD(lower = 0,upper = Inf, mu = muUi, Sigma = SigmaUi, dist="t",nu=nu)$varcov


          U0<-as.numeric(auxU1/auxU2)/auxcte
          U1<-(U0)*(MoMT$mean)
          U2<-(U0)*(MoMT$EYY)

          Auxtuy<-(matrix(y[i],2,1))

          uy<- Auxtuy*U0
          uy[2]<- U1

          uyy<-(Auxtuy%*%t(Auxtuy))

          AAx<-uyy[1,1]*U0
          ABx<-Auxtuy[1]%*%t(U1)
          BBx<-U2

          uyy[1,1]<- AAx
          uyy[1,2]<- ABx
          uyy[2,1]<- ABx
          uyy[2,2]<- BBx
   		}
      else{
        SigmaUi<- Sigma
        SigmaUiA<- Sigma*nu/(nu+2)
        auxupper<- -mu2[i]

        auxU1 <- stats::pt(auxupper/sqrt(SigmaUiA[2,2]),nu+2)

        auxU2 <- stats::pt(auxupper/sqrt(SigmaUi[2,2]),nu)

         if(auxU2 == 0) auxU2<- .Machine$double.xmin

        MomNT1<- MomTrunc::meanvarTMD(lower = c(-Inf,-Inf),upper = c(Inf,0),mu = c(mu1[i],mu2[i]),Sigma = SigmaUiA,dist="t",nu=nu+2)
        vary = MomTrunc::meanvarTMD(lower = c(-Inf,-Inf),upper = c(Inf,0),mu = c(mu1[i],mu2[i]), Sigma = Sigma,dist="t",nu=nu)$varcov

         U0<-as.numeric(auxU1/auxU2)
         U1<-auxU1/auxU2*MomNT1$mean
         U2<-auxU1/auxU2*MomNT1$EYY

         uy<-U1
        uyy<-U2

      }
     SigI <- solve(Sigma,tol=1e-30)
     aux41<- rbind(c(x[i,],0*w[i,]),c(x[i,]*0,w[i,]))
     aux5 <- aux41%*%betaC
     aux6<-  uyy-aux5%*%t(uy)- t(aux5%*%t(uy))+U0*aux5%*%t(aux5)
     suma1<- suma1+ aux6[1,1]-rhoa*(aux6[1,2]+aux6[2,1])+aux6[2,2]*rhoa^2
     suma11<-suma11+(aux6[1,2]+aux6[2,1])/2
     suma12<-suma12+aux6[2,2]

     suma2<- suma2+ U0*t(aux41)%*%SigI%*%(aux41)
     suma3<- suma3+ t(aux41)%*%SigI%*%uy
  }
  suma2aux<-(suma2+t(suma2))/2
  betaC<- solve(suma2aux, tol=1e-30)%*%suma3
  psi<- suma1/n
  rhoa<- suma11/suma12
  sigma2<- psi+rhoa^2
  sigma<- sqrt(sigma2)
  rho<- rhoa/sigma
  Sigma<- matrix(c(sigma2,rhoa,rhoa,1),2,2)

  beta<- betaC[1:p]
  gama<- betaC[(p+1):(p+q)]

 lkante1<- lkante
 ## estimating nu
 maxi<- stats::optimize(likeLt, c(3.001,150), tol = 1e-07, maximum = TRUE, y, x, w, cc, beta, gama, Sigma)
 nu<- maxi$maximum
 lkante<-maxi$objective

 criterio<-sqrt(abs(1-lkante1/lkante))

  if((rho >=1 - 1e-07)&&(rho <= -1 + 1e-07)) {
     stop("EM did not coverge")#criterio<-0.0000000000001
  }
 }

 if(verbose) cat("\n")
 AICc<- NULL
 AICcorr <- NULL
 BICc <- NULL
 desvios <- NULL

 out <- list(y=y, x=x, w=w, cc=cc, beta=beta, gamma=gama, rho=rho, sigma=sigma, nu=nu, sd=desvios, logL=lkante,AIC=AICc, AICc=AICcorr, BIC=BICc, family="T")
 class(out) <- "HeckmanEM"

 if(criteria){
  crit <- HeckmanEM.criteria(out)
  out$AIC <- crit$AIC
  out$AICc <- crit$AICc
  out$BIC <- crit$BIC
 }

 if(im){
  desvios <- HeckmanEM.infomat(out)
  out$sd <- desvios
 }

 if(round(nu,5) == 150) nu <- Inf

 return(out)
}



################################################################################
### Likelihood function  : Contaminated Normal
################################################################################

likeLcn <- function(y, x, w, cc, beta, gama, Sigma, nu1=0.1, nu2=0.1){
  sigma2  <- Sigma[1,1]
  sigma   <- sqrt(sigma2)
  rho     <- Sigma[1,2]/sigma
  med     <- (y - x%*%beta)/sigma
  mu.c    <- w%*%gama + rho*med
  sigma.c <- sqrt(1 - rho^2)
  aux0 <- dnorm(med*sqrt(nu2), log = T)- log(sigma)+0.5*log(nu2)
  aux01<- dnorm(med, log = T)- log(sigma)
  aux02<- nu1*exp(aux0)+(1-nu1)*exp(aux01)
  wp<- nu1*exp(aux0)/aux02
  aux1<- log(aux02) + log(wp*exp(pnorm(q = mu.c*sqrt(nu2)/sigma.c, log.p = FALSE))+(1-wp)*exp(pnorm(q = mu.c/sigma.c, log.p = FALSE)))
  aux2<- log(nu1*exp(pnorm(q = -w%*%gama*sqrt(nu2), log.p = FALSE))+(1-nu1)*exp(pnorm(q = -w%*%gama, log.p = FALSE)))
  fd      <- ifelse(cc==1,aux1 , aux2)
  return(sum(fd))
}


################################################
######contaminated normal density ##############
################################################
den.nc<-function(y, mu, sigma, nu1, nu2){
  dens<-nu1*dnorm(y,mu,sigma/sqrt(nu2))+(1-nu1)*dnorm(y,mu,sigma)
  return(dens)
}



###############################
# Cotaminated normal algorithm##
################################

EMcn.alg <- function(y, x, w, cc, error = 1e-6,nu,iter.max = 5000,criteria,im,verbose){
  start.time <- Sys.time()

  n      <- nrow(x)
  y      <- matrix(y, n, 1)
  p      <- ncol(x)
  q      <- ncol(w)

  ini <- init.values(y,x,w,cc)


  beta <- ini$beta
  gama <- ini$gamma
  rho <- ifelse(abs(ini$rho) >= 1-1e-07, 0.5, ini$rho)
  sigma <-  ifelse(ini$sigma <= 0, 1, ini$sigma)
  nu1    <- nu[1] ##
  nu2    <- nu[2] ##
  sigma2 <- sigma^2
  rhoa   <- sigma*rho
  Sigma  <- matrix(c(sigma2, rhoa, rhoa, 1), ncol = 2)

  betaC  <- c(beta,gama)

  criterio <- 1
  count    <- 0
  lista    <- vector("list", n)

  lkante   <- likeLcn(y, x, w, cc, beta, gama, Sigma, nu1, nu2) ###

  while((criterio > error) && (count <= iter.max)){
    count  <- count + 1
    if(verbose) cat("Iteration: ", count,"\r")
    #print(count)
    mu1    <- x%*%beta
    mu2    <- w%*%gama
    SigI   <- solve(Sigma, tol = 1e-30) # Inversa de Sigma

    suma1  <- 0
    suma11 <- 0
    suma12 <- 0
    suma2  <- matrix(0, p+q, p+q)
    suma3  <- matrix(0, p+q, 1)
    sumaep <- 0
    sumaE1 <-matrix(0, 2, 2)
    for (i in 1:n){
      #print(i)
      uy   <- epy<-  matrix(y[i], 2, 1)
      uyy  <- epyy<- matrix(0   , 2, 2)

      if(cc[i]==1){
        mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
        sigma12.1 <- 1 - rho^2
        omega.nu2 <- nu1*dnorm(y[i], mu1[i], sigma/sqrt(nu2))/den.nc(y[i],mu1[i], sigma, nu1, nu2)###

        if(sigma12.1 == 0) sigma12.1 <- 0.0001
        #' @importFrom MomTrunc meanvarTMD
        MomNT     <- MomTrunc::meanvarTMD(lower = 0, upper = Inf, mu = mu12.1, Sigma = sigma12.1, dist = "normal",n = 10^5)
        MomNT1    <- MomTrunc::meanvarTMD(lower = 0, upper = Inf, mu = mu12.1, Sigma = sigma12.1/nu2, dist = "normal",n = 10^5) ###

        ep0       <- pnorm(0, mu12.1, sqrt(sigma12.1/nu2),lower.tail = FALSE)
        ep1       <- pnorm(0, mu12.1, sqrt(sigma12.1),lower.tail = FALSE)

        epaux     <- omega.nu2*ep0+(1-omega.nu2)*ep1

        ep        <- omega.nu2*ep0/epaux

        wc1       <- (ep0*omega.nu2*MomNT1$mean+ep1*(1-omega.nu2)*MomNT$mean)/epaux ###
        wc2       <- (ep0*omega.nu2*MomNT1$EYY+ep1*(1-omega.nu2)*MomNT$EYY)/epaux ###

        w1c1      <- MomNT1$mean   ###
        w1c2      <- MomNT1$EYY  ###
        uy[2]     <- wc1 ###
        uyy       <- uy%*%t(uy) ###
        uyy[2,2]  <- wc2  ###

        epy[2]    <- w1c1
        epyy      <- epy%*%t(epy)
        epyy[2,2] <- w1c2
        epyy      <- epyy*ep
        epy       <- epy*ep
      }

      else{
        MomNT2    <- MomTrunc::meanvarTMD(lower = c(-Inf, -Inf), upper = c(Inf, 0), mu = c(mu1[i], mu2[i]), Sigma = Sigma, dist = "normal",n = 10^5)
        MomNT3    <- MomTrunc::meanvarTMD(lower = c(-Inf, -Inf), upper = c(Inf, 0), mu = c(mu1[i], mu2[i]), Sigma = Sigma/nu2, dist = "normal",n = 10^5)
        #' @importFrom mvtnorm pmvnorm
        auxp1      <- pmvnorm(lower=-Inf,upper=c(Inf,0), mean=c(mu1[i], mu2[i]), sigma =Sigma/nu2)[1]
        auxp2      <- pmvnorm(lower=-Inf,upper=c(Inf,0), mean=c(mu1[i], mu2[i]), sigma =Sigma)[1]
        auxp3      <- nu1*auxp1+(1-nu1)*auxp2

        Wc1       <- (auxp1*nu1*MomNT3$mean+auxp2*(1-nu1)*MomNT2$mean)/auxp3 ###
        Wc2       <- (auxp1*nu1*MomNT3$EYY+auxp2*(1-nu1)*MomNT2$EYY)/auxp3 ###

        W1c1      <- MomNT3$mean   ###
        W1c2      <- MomNT3$EYY  ###

        uy        <- Wc1
        uyy       <- Wc2


        ep        <- nu1*auxp1/auxp3
        epy       <- ep*W1c1
        epyy      <- ep*W1c2
      }

      aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))                          # X_{ic}
      aux5        <- aux41%*%betaC                                                        # mu_i = X_{ic} * betaC
      aux6.1      <- uyy - aux5%*%t(uy) - t(aux5%*%t(uy)) + aux5%*%t(aux5)                # Matrix Gamma_i
      aux6.2      <- epyy - aux5%*%t(epy) - t(aux5%*%t(epy)) + ep*aux5%*%t(aux5)
      sumaE1      <- sumaE1 + aux6.2
      aux6        <- aux6.1+(nu2-1)*aux6.2
      suma1       <- suma1  + aux6[1,1] - rhoa*(aux6[1,2] + aux6[2,1]) + aux6[2,2]*rhoa^2 # Cálculo de psi
      suma11      <- suma11 + (aux6[1,2] + aux6[2,1])/2                                   # Cálculo do numerador de rho*
      suma12      <- suma12 + aux6[2,2]                                                   # Cálculo do denominador de rho*
      suma2       <- suma2  + (1+(nu2-1)*ep)*t(aux41)%*%SigI%*%(aux41)                    # Cálculo do betaC (first sum)
      suma3       <- suma3  + t(aux41)%*%SigI%*%(uy+(nu2-1)*epy)                          # Cálculo do betaC (second sum)
      sumaep      <- sumaep + ep
      lista[[i]]        <- list(uy, uyy)
      names(lista[[i]]) <- c("E11", "E12")
    }

    suma2aux      <- (suma2 + t(suma2))/2 # Para forçar a diagnonal secundária ter os mesmos valores

    betaC         <- solve(suma2aux, tol = 1e-30)%*%suma3
    psi           <- suma1/n
    rhoa          <- suma11/suma12
    sigma2        <- psi + rhoa^2
    sigma         <- sqrt(sigma2)
    rho           <- rhoa/sigma
    Sigma         <- matrix(c(sigma2, rhoa, rhoa, 1), 2, 2)
    nu1           <- sumaep/n
    #print(nu1)
    nu2           <- min(1,2*n*nu1/sum(diag(solve(Sigma, tol = 1e-30)%*%sumaE1)))

    beta          <- betaC[1:p]
    gama          <- betaC[(p+1):(p+q)]

    lkante1       <- lkante
    lkante        <- likeLcn(y, x, w, cc, beta, gama, Sigma, nu1, nu2)
    criterio      <- sqrt(abs(1 - lkante1/lkante))

    end.time      <- Sys.time()
    time.inter    <- end.time - start.time

    #cat("Interaction: ", count, "Max. error: ", error,
    #   "Convergence: ", round(criterio, 6),
    #  "Time: ", round(time.inter,2),"\n")

    if(rho <= -0.999 | rho >= 0.999) stop("EM did not coverge") #criterio <- 1e-15
  }

  if(verbose) cat("\n")
  AICc<- NULL
  AICcorr <- NULL
  BICc <- NULL
  desvios <- NULL


  # theta_novo     <- matrix(c(beta, gama, sigma2, rho, nu1, nu2), ncol = 1)

  #npar           <- length(p+q+4)

  # Model comparison criteria
  #AICc           <- -2*lkante + 2*npar
  #AICcorr        <-    AICc   + ((2*npar*(npar + 1))/(n - npar - 1))
  #BICc           <- -2*lkante + log(n)*npar

  #desvios <- matrizMI(y, x, w, cc, beta, gama, rho, sigma, nu = c(nu1,nu2), type = "CN")

  out <- list(y=y, x=x, w=w, cc=cc, beta=beta, gamma=gama, rho=rho, sigma=sigma, nu=c(nu1,nu2), sd=desvios, logL=lkante,AIC=AICc, AICc=AICcorr, BIC=BICc, family="CN")
  class(out) <- "HeckmanEM"

  if(criteria){
    crit <- HeckmanEM.criteria(out)
    out$AIC <- crit$AIC
    out$AICc <- crit$AICc
    out$BIC <- crit$BIC
  }

  if(im){
    desvios <- HeckmanEM.infomat(out)
    out$sd <- desvios
  }

  #return(list(theta=theta_novo, beta=beta, gama=gama, rho=rho, sigma=sigma, nu1=nu1, nu2=nu2, sd=desvios, iter=count,
  #          logL=lkante, AIC=AICc, AICc=AICcorr, BIC=BICc, lista=lista))
  return(out)
}


