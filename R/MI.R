#' Standard error estimation for the Heckman Selection model by the Information Matrix
#'
#' `HeckmanEM.infomat()` estimates the standard errors for the parameters for the fitted Heckman selection model.
#'
#' @param obj An object of the class HeckmanEM.
#' @return The estimated standard errors for the parameters of the fitted model.
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
#'                  im = FALSE, criteria = TRUE)
#' im <- HeckmanEM.infomat(res)
#' }
#' @export
HeckmanEM.infomat <- function(obj){

  if(!inherits(obj,"HeckmanEM")) stop("Only \"HeckmanEM\" objects accepted!")
  if (obj$family != "Normal" && obj$family !="normal" && obj$family !="T" && obj$family !="t"  &&
      obj$family !="CN" && obj$family !="cn") stop("Family not recognized")

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

  n<-nrow(x)
  y<-matrix(y,n,1)
  p<-ncol(x)
  q<-ncol(w)

  sigma2 <- sigma^2
  rhoa<- sigma*rho
  Sigma<- matrix(c(sigma2, rhoa, rhoa,1), ncol = 2)

  betaC<-c(beta,gama)

  if (family=="Normal" || family=="normal"){

    mu1<- x%*%beta
    mu2<- w%*%gama
    suma1<-matrix(0,p+q+2,p+q+2)

    for (i in 1:n){

      uy <-matrix(y[i],2,1)
      uyy<-matrix(0,2,2)

      if(cc[i]==1){
        mu12.1<- mu2[i]+rho/sigma*(y[i]-mu1[i])
        sigma12.1<- 1-rho^2
        MomNT<- MomTrunc::meanvarTMD(lower = 0,upper = Inf, mu = mu12.1, Sigma = sigma12.1, dist="normal")

        uy[2]<-MomNT$mean
        uyy[2,2]<-MomNT$varcov
        uyy<- uyy+uy%*%t(uy)

      }
      else{

        MomNT1<- MomTrunc::meanvarTMD(lower = c(-Inf,-Inf),upper = c(Inf,0),mu = c(mu1[i],mu2[i]),Sigma = Sigma,dist="normal")
        uy <- MomNT1$mean
        uyy <- MomNT1$EYY

      }
      aux41<- rbind(c(x[i,],0*w[i,]),c(x[i,]*0,w[i,]))

      aux5 <- aux41%*%betaC
      aux2<-  uyy-aux5%*%t(uy)- t(aux5%*%t(uy))+aux5%*%t(aux5)

      Baa<-matrix(c(2*sigma, rho, rho ,0),ncol=2,nrow=2)
      Caa<-matrix(c(0, sigma,sigma,0),nrow=2,ncol=2)

      SigI <- chol2inv(chol(Sigma))
      derBetas<- 0.5*(t(aux41)%*%SigI%*%(uy))+0.5*t(t(uy)%*%SigI%*%(aux41))-t(aux41)%*%SigI%*%aux5
      dersigma<- -0.5*sum(diag(SigI%*%Baa))+0.5*sum(diag(aux2%*%SigI%*%Baa%*%SigI))
      derrho<- -0.5*sum(diag(SigI%*%Caa))+0.5*sum(diag(aux2%*%SigI%*%Caa%*%SigI))

      derlog<-  matrix(c(derBetas,dersigma,derrho),p+q+2,1)
      suma1<- suma1+derlog%*%t(derlog)
    }

  }

  if (family=="T" || family=="t"){


    ychap <- cbind(rep(1,n), rep(1,n))
    vary <- diag(2)

    mu1<- x%*%beta
    mu2<- w%*%gama


    suma1<-matrix(0,p+q+2,p+q+2)

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

        auxU1<- 1-stats::pt(auxupper/sqrt(SigmaUiA),nu1+2)

        auxU2<- 1-stats::pt(auxupper/sqrt(SigmaUi),nu1)

        MoMT<- MomTrunc::meanvarTMD(lower = 0,upper = Inf,mu = muUi, Sigma = SigmaUiA, dist="t",nu=nu1+2)


        vary <- matrix(0,2,2)
        vary[2,2]<-MomTrunc::meanvarTMD(lower = 0,upper = Inf, mu = muUi, Sigma = SigmaUi, dist="t",nu=nu)$varcov


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

        MomNT1<- MomTrunc::meanvarTMD(lower = c(-Inf,-Inf),upper = c(Inf,0),mu = c(mu1[i],mu2[i]),Sigma = SigmaUiA,dist="t",nu=nu+2)
        vary = MomTrunc::meanvarTMD(lower = c(-Inf,-Inf),upper = c(Inf,0),mu = c(mu1[i],mu2[i]), Sigma = Sigma,dist="t",nu=nu)$varcov

        U0<-as.numeric(auxU1/auxU2)
        U1<-auxU1/auxU2*MomNT1$mean
        U2<-auxU1/auxU2*MomNT1$EYY

        uy<-U1
        uyy<-U2

      }

      aux41<- rbind(c(x[i,],0*w[i,]),c(x[i,]*0,w[i,]))

      aux5 <- aux41%*%betaC
      aux2<-  uyy-aux5%*%t(uy)- t(aux5%*%t(uy))+U0*aux5%*%t(aux5)

      Baa<-matrix(c(2*sigma, rho, rho ,0),ncol=2,nrow=2)
      Caa<-matrix(c(0, sigma,sigma,0),nrow=2,ncol=2)

      SigI <- chol2inv(chol(Sigma))
      derBetas<- 0.5*(t(aux41)%*%SigI%*%(uy))+0.5*t(t(uy)%*%SigI%*%(aux41))-U0*t(aux41)%*%SigI%*%aux5
      dersigma<- -0.5*sum(diag(SigI%*%Baa))+0.5*sum(diag(aux2%*%SigI%*%Baa%*%SigI))
      derrho<- -0.5*sum(diag(SigI%*%Caa))+0.5*sum(diag(aux2%*%SigI%*%Caa%*%SigI))

      derlog<-  matrix(c(derBetas,dersigma,derrho),p+q+2,1)
      suma1<- suma1+derlog%*%t(derlog)
    }

  }

  if (family =="CN" || family == "cn"){
    mu1   <- x%*%beta
    mu2   <- w%*%gama
    nu1<-nu[1]
    nu2<-nu[2]
    suma1 <- matrix(0, p+q+4, p+q+4)

    for (i in 1:n){
      uy   <- epy<-  matrix(y[i], 2, 1)
      uyy  <- epyy<- matrix(0   , 2, 2)

      if(cc[i]==1){
        mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
        sigma12.1 <- 1 - rho^2
        omegaxnu1 = dnorm(y[i], mu1[i], sigma/sqrt(nu2))/den.nc(y[i],mu1[i], sigma, nu1, nu2)
        omega.nu2 <- nu1*omegaxnu1###

        if(sigma12.1 == 0) sigma12.1 <- 0.0001

        MomNT     <- MomTrunc::meanvarTMD(lower = 0, upper = Inf, mu = mu12.1, Sigma = sigma12.1, dist = "normal",n = 10^5)
        MomNT1    <- MomTrunc::meanvarTMD(lower = 0, upper = Inf, mu = mu12.1, Sigma = sigma12.1/nu2, dist = "normal",n = 10^5) ###

        ep0       <- pnorm(0, mu12.1, sqrt(sigma12.1/nu2),lower.tail = FALSE)
        ep1       <- pnorm(0, mu12.1, sqrt(sigma12.1),lower.tail = FALSE)

        epaux     <- omega.nu2*ep0+(1-omega.nu2)*ep1

        ep        <- omega.nu2*ep0/epaux
        epovernu1 = omegaxnu1*ep0/epaux
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

        auxp1      <- pmvnorm(lower=-Inf,upper=c(Inf,0), mean=c(mu1[i], mu2[i]), sigma =Sigma/nu2)[1]
        auxp2      <- pmvnorm(lower=-Inf,upper=c(Inf,0), mean=c(mu1[i], mu2[i]), sigma =Sigma)[1]
        auxp3      <- nu1*auxp1+(1-nu1)*auxp2

        Wc1       <- (auxp1*nu1*MomNT3$mean+auxp2*(1-nu1)*MomNT2$mean)/auxp3 ###
        Wc2       <- (auxp1*nu1*MomNT3$EYY+auxp2*(1-nu1)*MomNT2$EYY)/auxp3 ###

        W1c1      <- MomNT3$mean   ###
        W1c2      <- MomNT3$EYY  ###

        uy        <- Wc1
        uyy       <- Wc2

        epovernu1 =  auxp1/auxp3
        ep        <- nu1*auxp1/auxp3
        epy       <- ep*W1c1
        epyy      <- ep*W1c2
      }

      aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))

      aux5        <- aux41%*%betaC
      aux2.1      <- uyy - aux5%*%t(uy) - t(aux5%*%t(uy)) + aux5%*%t(aux5)
      aux2.2      <- epyy - aux5%*%t(epy) - t(aux5%*%t(epy)) + ep*aux5%*%t(aux5)
      aux2        <- aux2.1+(nu2-1)*aux2.2

      Baa         <- matrix(c(2*sigma, rho  , rho  , 0), ncol = 2, nrow = 2)
      Caa         <- matrix(c(0      , sigma, sigma, 0), nrow = 2, ncol = 2)

      SigI        <- chol2inv(chol(Sigma))
      derBetas1    <-  0.5*(t(aux41)%*%SigI%*%(uy)) + 0.5*t(t(uy)%*%SigI%*%(aux41)) - t(aux41)%*%SigI%*%aux5
      derBetas2    <-  (nu2-1)*(0.5*(t(aux41)%*%SigI%*%(epy)) + 0.5*t(t(epy)%*%SigI%*%(aux41)) - ep*t(aux41)%*%SigI%*%aux5)
      derBetas   <-derBetas1+derBetas2
      dersigma    <- -0.5*sum(diag(SigI%*%Baa))    + 0.5*sum(diag(aux2%*%SigI%*%Baa%*%SigI))
      derrho      <- -0.5*sum(diag(SigI%*%Caa))    + 0.5*sum(diag(aux2%*%SigI%*%Caa%*%SigI))
      dernu1      <-  epovernu1 +(1-ep)/(1-nu1)
      dernu2      <-  ep/nu2-0.5* sum(diag(SigI%*%aux2.2))
      derlog      <- matrix(c(derBetas, dersigma, derrho,dernu1,dernu2), p+q+4, 1)

      suma1       <- suma1 + derlog%*%t(derlog)
    }
  }


  if (family =="CN" || family == "cn"){
    s = nrow(suma1)
    matvar1 = chol2inv(chol(suma1[1:(s-1),1:(s-1)]))
    if (det(suma1)==0){
      desvio = sqrt(diag(matvar1))
      names(desvio) <- c(paste0("beta",0:(length(beta)-1)), paste0("gamma",0:(length(gama)-1)), "sigma", "rho","nu1")
      warning("It was not possible to compute the standard error for nu2, ill-conditioned inverse matrix")
    }else{
      matvar = chol2inv(chol(suma1))
      if(matvar[s,s]>=1){
        desvio = sqrt(diag(matvar1))
        names(desvio) <- c(paste0("beta",0:(length(beta)-1)), paste0("gamma",0:(length(gama)-1)), "sigma", "rho","nu1")
        warning("It was not possible to compute the standard error for nu2, ill-conditioned inverse matrix")
      }else{
        desvio = sqrt(diag(matvar))
        names(desvio) <- c(paste0("beta",0:(length(beta)-1)), paste0("gamma",0:(length(gama)-1)), "sigma", "rho","nu1","nu2")
      }

    }
  }else{
    desvio<- sqrt(diag(chol2inv(chol(suma1))))
    names(desvio) <- c(paste0("beta",0:(length(beta)-1)), paste0("gamma",0:(length(gama)-1)), "sigma", "rho")
  }


  return(desvio)
}
