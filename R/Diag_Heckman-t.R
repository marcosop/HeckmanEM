### ==========================================================
### Influence diagnostics in Heckman selection-$t$ model
### Authors: Marcos S. Oliveira and Victor H. Lachos
### Last version: 02/12/2023
### ==========================================================

## -------------------------------
## Q FUNCTION
## -------------------------------

Q.function <- function(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, type = "Normal"){

  n        <- nrow(x)
  y        <- matrix(y, n, 1)
  p        <- ncol(x)
  q        <- ncol(w)

  sigma2   <- sigma^2
  rhoa     <- sigma*rho
  Sigma    <- matrix(c(sigma2, rhoa, rhoa, 1), ncol = 2)
  betaC    <- c(beta,gama)

  if (type=="Normal"){
    mu1   <- x%*%beta
    mu2   <- w%*%gama
    suma  <- matrix(0, 1, 1)

    for (i in 1:n){
      uy  <- matrix(y[i],2,1)
      uyy <- matrix(0,2,2)

      if(cc[i]==1){
        mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
        sigma12.1 <- 1 - rho^2
        MomNT     <- meanvarTMD(0, Inf, mu12.1, sigma12.1, dist = "normal")

        uy[2]     <- MomNT$mean
        uyy[2,2]  <- MomNT$varcov
        uyy       <- uyy + uy%*%t(uy)
      }
      else{
        MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma, dist = "normal")
        uy        <- MomNT1$mean
        uyy       <- MomNT1$EYY
      }

      aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
      aux5        <- aux41%*%betaC                                      # Xic * Bc
      aux2        <- uyy - uy%*%t(aux5) - aux5%*%t(uy) + aux5%*%t(aux5) # Gamma_i

      SigI        <- solve(Sigma)

      suma        <- suma - 0.5*log(det(Sigma)) - 0.5*sum(diag(aux2%*%SigI))
    }
  }

  if (type=="T"){
    ychap <- cbind(rep(1,n), rep(1,n))
    vary  <- diag(2)

    mu1   <- x%*%beta
    mu2   <- w%*%gama

    suma  <- matrix(0, 1, 1)

    for (i in 1:n){
      uy  <- matrix(y[i],2,1)
      uyy <- matrix(0,2,2)

      if(cc[i]==1){
        PsiA      <- Sigma*nu/(nu+2)
        nu1       <- (nu + 1)                            # X
        mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
        sigma12.1 <- 1 - rho^2
        ScA       <- nu/(nu + 2)*sigma12.1

        Qy1       <- (y[i] - mu1[i])^2/sigma2            # X
        Qy2       <- (y[i] - mu1[i])^2/PsiA[1, 1]        # X
        auxcte    <- as.numeric((nu + Qy1)/(nu + 1))
        auxcte1   <- as.numeric((nu + 2 + Qy2)/(nu + 3))

        Sc22      <- auxcte*sigma12.1

        muUi      <- mu12.1
        SigmaUi   <- Sc22

        SigmaUiA  <- auxcte1*ScA
        auxupper  <- -muUi

        auxU1     <- 1 - pt(auxupper/sqrt(SigmaUiA), nu1 + 2)
        auxU2     <- 1 - pt(auxupper/sqrt(SigmaUi), nu1)

        MoMT      <- meanvarTMD(0, Inf, muUi, SigmaUiA, dist = "t", nu = nu1 + 2)

        vary      <- matrix(0, 2, 2)
        vary[2,2] <- meanvarTMD(0, Inf, muUi, SigmaUi , dist = "t", nu = nu)$varcov

        U0        <- as.numeric(auxU1/auxU2)/auxcte
        U1        <- (U0)*(MoMT$mean)
        U2        <- (U0)*(MoMT$EYY)

        Auxtuy    <- (matrix(y[i], 2, 1))

        uy        <- Auxtuy*U0
        uy[2]     <- U1

        uyy       <- (Auxtuy%*%t(Auxtuy))

        AAx       <- uyy[1, 1]*U0
        ABx       <- Auxtuy[1]%*%t(U1)
        BBx       <- U2

        uyy[1,1]  <- AAx
        uyy[1,2]  <- ABx
        uyy[2,1]  <- ABx
        uyy[2,2]  <- BBx
      }
      else{
        SigmaUi   <- Sigma
        SigmaUiA  <- Sigma*nu/(nu+2)
        auxupper  <- -mu2[i]

        auxU1     <- pt(auxupper/sqrt(SigmaUiA[2,2]), nu+2)
        auxU2     <- pt(auxupper/sqrt(SigmaUi[2,2]) , nu)

        MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), SigmaUiA, dist = "t", nu = nu+2)
        vary      <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma   , dist = "t", nu = nu)$varcov

        U0        <- as.numeric(auxU1/auxU2)
        U1        <- auxU1/auxU2*MomNT1$mean
        U2        <- auxU1/auxU2*MomNT1$EYY

        uy        <- U1
        uyy       <- U2
      }

      aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))           # Xic
      aux5        <- aux41%*%betaC                                         # Xic * Bc
      aux2        <- uyy - uy%*%t(aux5) - aux5%*%t(uy) + U0*aux5%*%t(aux5) # Gamma_i atualizado

      SigI        <- solve(Sigma)

      suma        <- suma - 0.5*log(det(Sigma)) - 0.5*sum(diag(aux2%*%SigI))
    }
  }

  return(suma)
}

### -------------------------------
### HESSIAN Q
### -------------------------------

HessianaQ <- function(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, type = "Normal"){

  n        <- nrow(x)
  y        <- matrix(y, n, 1)
  p        <- ncol(x)
  q        <- ncol(w)

  sigma2   <- sigma^2
  rhoa     <- sigma*rho
  Sigma    <- matrix(c(sigma2, rhoa, rhoa, 1), ncol = 2)
  betaC    <- c(beta,gama)

  d.Ssigma2     <- matrix(c(1, rho/(2*sigma), rho/(2*sigma), 0), 2, 2)
  d.Srho        <- matrix(c(0, sigma, sigma, 0), 2, 2)

  d2.Ssigma2    <- matrix(c(0, -rho/(4*sigma^3), -rho/(4*sigma^3), 0), 2, 2)
  d2.Ssigma2rho <- matrix(c(0, 1/(2*sigma), 1/(2*sigma), 0), 2, 2)
  d2.Srho       <- matrix(0, 2, 2)

  theta    <- matrix(c(beta, gama, sigma2, rho), ncol = 1)

  if (type=="Normal"){
    mu1     <- x%*%beta
    mu2     <- w%*%gama

    suma1   <- matrix(0, p+q, p+q) # betacbetac
    suma2   <- matrix(0, p+q, 1)   # betacsigma2
    suma3   <- matrix(0, p+q, 1)   # betacrho
    suma4   <- 0                   # sigma2sigma2
    suma5   <- 0                   # sigma2rho
    suma6   <- 0                   # rhorho

    sbetac  <- matrix(0, n, p+q)
    ssigma2 <- matrix(0, n, 1)
    srho    <- matrix(0, n, 1)

    for (i in 1:n){
      uy  <- matrix(y[i],2,1)
      uyy <- matrix(0,2,2)

      if(cc[i]==1){
        mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
        sigma12.1 <- 1 - rho^2
        MomNT     <- meanvarTMD(0, Inf, mu12.1, sigma12.1, dist = "normal")

        uy[2]     <- MomNT$mean
        uyy[2,2]  <- MomNT$varcov
        uyy       <- uyy + uy%*%t(uy)
      }
      else{
        MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma, dist = "normal")
        uy        <- MomNT1$mean
        uyy       <- MomNT1$EYY
      }

      aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
      aux5        <- aux41%*%betaC                                      # Xic * Bc
      aux2        <- uyy - uy%*%t(aux5) - aux5%*%t(uy) + aux5%*%t(aux5) # Gamma_i

      SigI        <- solve(Sigma)

      suma1       <- suma1 - t(aux41)%*%SigI%*%aux41
      suma2       <- suma2 - t(aux41)%*%SigI%*%d.Ssigma2%*%SigI%*%uy + t(aux41)%*%SigI%*%d.Ssigma2%*%SigI%*%aux5
      suma3       <- suma3 - t(aux41)%*%SigI%*%d.Srho%*%SigI%*%uy + t(aux41)%*%SigI%*%d.Srho%*%SigI%*%aux5

      suma4       <- suma4 + 0.5*sum(diag(SigI%*%d.Ssigma2%*%SigI%*%d.Ssigma2 - SigI%*%d2.Ssigma2))
      + 0.5*sum(diag(aux2%*%(SigI%*%d2.Ssigma2%*%SigI - 2*SigI%*%d.Ssigma2%*%SigI%*%d.Ssigma2%*%SigI)))
      suma5       <- suma5 + 0.5*sum(diag(SigI%*%d.Srho%*%SigI%*%d.Ssigma2 - SigI%*%d2.Ssigma2rho))
      + 0.5*sum(diag(aux2%*%(SigI%*%d2.Ssigma2rho%*%SigI - SigI%*%d.Srho%*%SigI%*%d.Ssigma2%*%SigI - SigI%*%d.Ssigma2%*%SigI%*%d.Srho%*%SigI)))
      suma6       <- suma6 + 0.5*sum(diag(SigI%*%d.Srho%*%SigI%*%d.Srho - SigI%*%d2.Srho))
      + 0.5*sum(diag(aux2%*%(SigI%*%d2.Srho%*%SigI - 2*SigI%*%d.Srho%*%SigI%*%d.Srho%*%SigI)))

      sbetac[i,]  <- t(aux41)%*%SigI%*%uy - t(aux41)%*%SigI%*%aux5
      ssigma2[i]  <- -0.5*sum(diag(SigI%*%d.Ssigma2)) + 0.5*sum(diag(aux2%*%SigI%*%d.Ssigma2%*%SigI))
      srho[i]     <- -0.5*sum(diag(SigI%*%d.Srho))    + 0.5*sum(diag(aux2%*%SigI%*%d.Srho%*%SigI))
    }

    suma1aux         <- (suma1 + t(suma1))/2

    derbetacbetac    <- suma1aux
    derbetacsigma2   <- suma2
    derbetacrho      <- suma3
    dersigma2sigma2  <- suma4
    dersigma2rho     <- suma5
    derrhorho        <- suma6

    MatrizQ                  <- matrix(0, nrow = (p+q+2), ncol = (p+q+2))
    MatrizQ[1:(p+q),1:(p+q)] <- derbetacbetac
    MatrizQ[p+q+1,1:(p+q)]   <- t(derbetacsigma2)
    MatrizQ[1:(p+q),p+q+1]   <- derbetacsigma2
    MatrizQ[p+q+2,1:(p+q)]   <- t(derbetacrho)
    MatrizQ[1:(p+q),p+q+2]   <- derbetacrho
    MatrizQ[p+q+1,p+q+1]     <- dersigma2sigma2
    MatrizQ[p+q+1,p+q+2]     <- dersigma2rho
    MatrizQ[p+q+2,p+q+1]     <- dersigma2rho
    MatrizQ[p+q+2,p+q+2]     <- derrhorho

    dQtheta                  <- matrix(0, p+q+2, n)
    thetai                   <- matrix(0, p+q+2, n)

    HesI                     <- solve(-MatrizQ) # Inversa da Hessiana x (-1)

    for(i in 1:n){
      dQtheta[1:(p+q),i]     <- apply(sbetac[-i,], 2, sum)
      dQtheta[p+q+1,i]       <- sum(ssigma2[-i])
      dQtheta[p+q+2,i]       <- sum(srho[-i])
      thetai[,i]             <- theta[1:(p+q+2)] + HesI%*%dQtheta[,i]
    }
  }

  if (type=="T"){
    ychap <- cbind(rep(1,n), rep(1,n))
    vary  <- diag(2)

    mu1   <- x%*%beta
    mu2   <- w%*%gama

    suma1   <- matrix(0, p+q, p+q) # betacbetac
    suma2   <- matrix(0, p+q, 1)   # betacsigma2
    suma3   <- matrix(0, p+q, 1)   # betacrho
    suma4   <- 0                   # sigma2sigma2
    suma5   <- 0                   # sigma2rho
    suma6   <- 0                   # rhorho

    sbetac  <- matrix(0, n, p+q)
    ssigma2 <- matrix(0, n, 1)
    srho    <- matrix(0, n, 1)

    for (i in 1:n){
      uy  <- matrix(y[i],2,1)
      uyy <- matrix(0,2,2)

      if(cc[i]==1){
        PsiA      <- Sigma*nu/(nu+2)
        nu1       <- (nu + 1)                            # X
        mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
        sigma12.1 <- 1 - rho^2
        ScA       <- nu/(nu + 2)*sigma12.1

        Qy1       <- (y[i] - mu1[i])^2/sigma2            # X
        Qy2       <- (y[i] - mu1[i])^2/PsiA[1, 1]        # X
        auxcte    <- as.numeric((nu + Qy1)/(nu + 1))
        auxcte1   <- as.numeric((nu + 2 + Qy2)/(nu + 3))

        Sc22      <- auxcte*sigma12.1

        muUi      <- mu12.1
        SigmaUi   <- Sc22

        SigmaUiA  <- auxcte1*ScA
        auxupper  <- -muUi

        auxU1     <- 1 - pt(auxupper/sqrt(SigmaUiA), nu1 + 2)
        auxU2     <- 1 - pt(auxupper/sqrt(SigmaUi), nu1)

        MoMT      <- meanvarTMD(0, Inf, muUi, SigmaUiA, dist = "t", nu = nu1 + 2)

        vary      <- matrix(0, 2, 2)
        vary[2,2] <- meanvarTMD(0, Inf, muUi, SigmaUi , dist = "t", nu = nu)$varcov

        U0        <- as.numeric(auxU1/auxU2)/auxcte
        U1        <- (U0)*(MoMT$mean)
        U2        <- (U0)*(MoMT$EYY)

        Auxtuy    <- (matrix(y[i], 2, 1))

        uy        <- Auxtuy*U0
        uy[2]     <- U1

        uyy       <- (Auxtuy%*%t(Auxtuy))

        AAx       <- uyy[1, 1]*U0
        ABx       <- Auxtuy[1]%*%t(U1)
        BBx       <- U2

        uyy[1,1]  <- AAx
        uyy[1,2]  <- ABx
        uyy[2,1]  <- ABx
        uyy[2,2]  <- BBx
      }
      else{
        SigmaUi   <- Sigma
        SigmaUiA  <- Sigma*nu/(nu+2)
        auxupper  <- -mu2[i]

        auxU1     <- pt(auxupper/sqrt(SigmaUiA[2,2]), nu+2)
        auxU2     <- pt(auxupper/sqrt(SigmaUi[2,2]) , nu)

        MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), SigmaUiA, dist = "t", nu = nu+2)
        vary      <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma   , dist = "t", nu = nu)$varcov

        U0        <- as.numeric(auxU1/auxU2)
        U1        <- auxU1/auxU2*MomNT1$mean
        U2        <- auxU1/auxU2*MomNT1$EYY

        uy        <- U1
        uyy       <- U2
      }

      aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))           # Xic
      aux5        <- aux41%*%betaC                                         # Xic * Bc
      aux2        <- uyy - uy%*%t(aux5) - aux5%*%t(uy) + U0*aux5%*%t(aux5) # Gamma_i atualizado

      SigI        <- solve(Sigma)

      suma1       <- suma1 - U0*t(aux41)%*%SigI%*%aux41
      suma2       <- suma2 - t(aux41)%*%SigI%*%d.Ssigma2%*%SigI%*%uy + U0*t(aux41)%*%SigI%*%d.Ssigma2%*%SigI%*%aux5
      suma3       <- suma3 - t(aux41)%*%SigI%*%d.Srho%*%SigI%*%uy    + U0*t(aux41)%*%SigI%*%d.Srho%*%SigI%*%aux5

      suma4       <- suma4 + 0.5*sum(diag(SigI%*%d.Ssigma2%*%SigI%*%d.Ssigma2 - SigI%*%d2.Ssigma2))
      + 0.5*sum(diag(aux2%*%(SigI%*%d2.Ssigma2%*%SigI - 2*SigI%*%d.Ssigma2%*%SigI%*%d.Ssigma2%*%SigI)))
      suma5       <- suma5 + 0.5*sum(diag(SigI%*%d.Srho%*%SigI%*%d.Ssigma2 - SigI%*%d2.Ssigma2rho))
      + 0.5*sum(diag(aux2%*%(SigI%*%d2.Ssigma2rho%*%SigI - SigI%*%d.Srho%*%SigI%*%d.Ssigma2%*%SigI - SigI%*%d.Ssigma2%*%SigI%*%d.Srho%*%SigI)))
      suma6       <- suma6 + 0.5*sum(diag(SigI%*%d.Srho%*%SigI%*%d.Srho - SigI%*%d2.Srho))
      + 0.5*sum(diag(aux2%*%(SigI%*%d2.Srho%*%SigI - 2*SigI%*%d.Srho%*%SigI%*%d.Srho%*%SigI)))

      sbetac[i,]  <- t(aux41)%*%SigI%*%uy - U0*t(aux41)%*%SigI%*%aux5
      ssigma2[i]  <- -0.5*sum(diag(SigI%*%d.Ssigma2)) + 0.5*sum(diag(aux2%*%SigI%*%d.Ssigma2%*%SigI))
      srho[i]     <- -0.5*sum(diag(SigI%*%d.Srho))    + 0.5*sum(diag(aux2%*%SigI%*%d.Srho%*%SigI))
    }

    suma1aux         <- (suma1 + t(suma1))/2

    derbetacbetac    <- suma1aux
    derbetacsigma2   <- suma2
    derbetacrho      <- suma3
    dersigma2sigma2  <- suma4
    dersigma2rho     <- suma5
    derrhorho        <- suma6

    MatrizQ                  <- matrix(0, nrow = (p+q+2), ncol = (p+q+2))
    MatrizQ[1:(p+q),1:(p+q)] <- derbetacbetac
    MatrizQ[p+q+1,1:(p+q)]   <- t(derbetacsigma2)
    MatrizQ[1:(p+q),p+q+1]   <- derbetacsigma2
    MatrizQ[p+q+2,1:(p+q)]   <- t(derbetacrho)
    MatrizQ[1:(p+q),p+q+2]   <- derbetacrho
    MatrizQ[p+q+1,p+q+1]     <- dersigma2sigma2
    MatrizQ[p+q+1,p+q+2]     <- dersigma2rho
    MatrizQ[p+q+2,p+q+1]     <- dersigma2rho
    MatrizQ[p+q+2,p+q+2]     <- derrhorho

    dQtheta                  <- matrix(0, p+q+2, n)
    thetai                   <- matrix(0, p+q+2, n)

    HesI                     <- solve(-MatrizQ) # Inversa da Hessiana x (-1)

    for(i in 1:n){
      dQtheta[1:(p+q),i]     <- apply(sbetac[-i,], 2, sum)
      dQtheta[p+q+1,i]       <- sum(ssigma2[-i])
      dQtheta[p+q+2,i]       <- sum(srho[-i])
      thetai[,i]             <- theta[1:(p+q+2)] + HesI%*%dQtheta[,i]
    }
  }

  obj.out <- list(MatrizQ = MatrizQ, dQtheta = dQtheta, thetai = thetai)

  return(obj.out)
}

### -------------------------------
### CASE-DELETION - User
### -------------------------------

#' Case deletion analysis for Heckman selection model
#'
#' This function performs case deletion analysis based on a HeckmanEM object (not available for the
#' contaminated normal model).
#'
#' @param object A HeckmanEM object.
#'
#' @return A list of class \code{HeckmanEM.deletion} with a vector GD of dimension \eqn{n} (see details), and a benchmark value.
#'
#' @details This function uses the case deletion approach to study
#' the impact of deleting one or more observations from the dataset
#' on the parameters estimates, using the ideas of Cook (1977) and Zhu et.al. (2001).
#' The GD vector contains the generalized Cook distances
#' \deqn{\textrm{GD}^1_i = \dot{Q}_{[i]}(\widehat{\boldsymbol{\theta}} \mid \widehat{\boldsymbol{\theta}})^{\top} \left\{-\ddot{Q}(\widehat{\boldsymbol{\theta}} \mid \widehat{\boldsymbol{\theta}})\right\}^{-1}\dot{Q}_{[i]}(\widehat{\boldsymbol{\theta}} \mid \widehat{\boldsymbol{\theta}}),}
#' where \eqn{\dot{Q}_{[i]}(\widehat{\boldsymbol{\theta}}\mid \widehat{\boldsymbol{\theta}})} is the gradient vector after dropping the \eqn{i}th observation, and
#' \eqn{\ddot{Q}(\widehat{\boldsymbol{\theta}} \mid \widehat{\boldsymbol{\theta}})} is the Hessian matrix. The benchmark was adapted using the suggestion of Barros et al. (2010). We use \eqn{(2 \times \textrm{npar})/n} as the benchmark for the \eqn{\textrm{GD}_i}, with \eqn{\textrm{npar}} representing the number of estimated model parameters.
#'
#' @export CaseDeletion
#'
#'@examples
#' n    <- 100
#' nu   <- 3
#' cens <- 0.25
#'
#' set.seed(13)
#' w <- cbind(1, runif(n, -1, 1), rnorm(n))
#' x <- cbind(w[,1:2])
#' c <- qt(cens, df = nu)
#'
#' sigma2   <- 1
#' beta     <- c(1, 0.5)
#' gamma    <- c(1, 0.3, -.5)
#' gamma[1] <- -c * sqrt(sigma2)
#'
#' datas <- rHeckman(x, w, beta, gamma, sigma2, rho = 0.6, nu, family = "T")
#' y     <- datas$y
#' cc    <- datas$cc
#'\donttest{
#' heckmodel <- HeckmanEM(y, x, w, cc, family = "Normal", iter.max = 50)
#'
#' global <- CaseDeletion(heckmodel)
#' plot(global)
#'}
#' @references M. Barros, M. Galea, M. González, V. Leiva, Influence diagnostics in the Tobit censored response model, Statistical Methods & Applications 19 (2010) 379–397.
#'
#' R. D. Cook, Detection of influential observation in linear regression, Technometrics 19 (1977) 15–18.
#'
#' H. Zhu, S. Lee, B. Wei, J. Zhou, Case-deletion measures for models with incomplete data, Biometrika 88 (2001) 727–737.
#'

CaseDeletion <- function(object){

  if(is(object,"HeckmanEM")){
    if(object$family == "CN" || object$family == "cn"){
      stop("case deletion analysis is not available for the contaminated normal model")
    }
    CaseDele(y = object$y,x = object$x,w = object$w,cc = object$cc,
             beta=object$beta, gama=object$gamma,rho=object$rho,
             sigma=object$sigma,nu = object$nu,
             family = object$family)

  }else{
    stop("The object must be of the `HeckmanEM` class.")
  }
}


CaseDele <- function(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, family = "Normal"){
  #start.time <- Sys.time()

  HQ   <- HessianaQ(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, type = family)

  HesI <- solve(-HQ$MatrizQ) # Inversa da Hessiana x (-1)

  n    <- nrow(x)
  p    <- ncol(x)
  q    <- ncol(w)

  GD       <- matrix(0, n, 1)

  if(family=="Normal"){
    beta   <- HQ$thetai[1:p,]
    gama   <- HQ$thetai[(p+1):(p+q),]
    sigma2 <- HQ$thetai[p+q+1,]
    rho    <- HQ$thetai[p+q+2,]
    sigma  <- sqrt(HQ$thetai[p+q+1,])

    for (i in 1:n){
      GD[i]  <- abs(t(HQ$dQtheta[,i])%*%HesI%*%HQ$dQtheta[,i])

      # end.time      <- Sys.time()
      # time.inter    <- end.time - start.time
      #
      # cat("% completed: ", round((i/n)*100,2), "Time: ", round(time.inter,2),"\n")
    }

    ### Benchmark for GD by Massuia et al. (2015)
    benchmark_GD  <- 2*(p+q+2)/n
  }

  if(family=="T"){
    beta   <- HQ$thetai[1:p,]
    gama   <- HQ$thetai[(p+1):(p+q),]
    sigma2 <- HQ$thetai[p+q+1,]
    rho    <- HQ$thetai[p+q+2,]
    sigma  <- sqrt(HQ$thetai[p+q+1,])

    for (i in 1:n){
      GD[i]  <- abs(t(HQ$dQtheta[,i])%*%HesI%*%HQ$dQtheta[,i])

      # end.time      <- Sys.time()
      # time.inter    <- end.time - start.time

      #cat("% completed: ", round((i/n)*100,2), "Time: ", round(time.inter,2),"\n")

      ### Benchmark for GD by Massuia et al. (2015)
      benchmark_GD  <- 2*(p+q+3)/n
    }
  }

  obj.out = list(GD = GD,benchmark = benchmark_GD)

  class(obj.out) = "HeckmanEM.deletion"

  return(obj.out)
}

### -------------------------------
### INFLUENCE DIAGNOSTICS - User
### -------------------------------

#' Influence Analysis for the Heckman Selection model
#'
#' @description
#' This function conducts influence analysis for a given `HeckmanEM` object.
#' The influence analysis can be conducted using several types of perturbations (not available for the
#' contaminated Normal model).
#'
#' @param object A `HeckmanEM` object to perform the analysis on.
#' @param type A character string indicating the type of perturbation to perform.
#' The types can be one of "case-weight","scale","response" and"exploratory".
#' @param colx Optional integer specifying the position of the column in the
#' object's matrix \code{x} that will undergo perturbation. Only required when type is
#' "exploratory".
#' @param k A positive real constant to be used in the benchmark calculation: \eqn{M_0 + k\times \mathrm{sd}(M_0)}. Default is 3.5.
#'
#' @return Returns a list of class \code{HeckmanEM.influence} with the following elements:
#' \item{M0}{A vector of length \eqn{n} with the aggregated contribution of all eigenvectors of the matrix associated with the normal curvature.}
#' \item{benchmark}{\eqn{M_0 + k\times \mathrm{sd}(M_0)}}
#' \item{influent}{A vector with the influential observations' positions.}
#' \item{type}{The perturbation type.}
#'
#' @examples
#' n    <- 100
#' nu   <- 3
#' cens <- 0.25
#'
#' set.seed(13)
#' w <- cbind(1, runif(n, -1, 1), rnorm(n))
#' x <- cbind(w[,1:2])
#' c <- qt(cens, df = nu)
#'
#' sigma2   <- 1
#' beta     <- c(1, 0.5)
#' gamma    <- c(1, 0.3, -.5)
#' gamma[1] <- -c * sqrt(sigma2)
#'
#' datas <- rHeckman(x, w, beta, gamma, sigma2, rho = 0.6, nu, family = "T")
#' y     <- datas$y
#' cc    <- datas$cc
#'\donttest{
#' heckmodel <- HeckmanEM(y, x, w, cc, family = "Normal", iter.max = 50)
#'
#' global <- CaseDeletion(heckmodel)
#' plot(global)
#'
#' local_case <- Influence(heckmodel, type = "case-weight")
#' local_case$influent # influential values here!
#' plot(local_case)
#'
#' local_scale <- Influence(heckmodel, type = "scale")
#' local_scale$influent # influential values here!
#' plot(local_scale)
#'
#' local_response <- Influence(heckmodel, type = "response")
#' local_response$influent # influential values here!
#' plot(local_response)
#'
#' local_explore <- Influence(heckmodel, type = "exploratory", colx = 2)
#' local_explore$influent # influential values here!
#' plot(local_explore)
#'
#'}
#' @seealso
#' \code{\link{HeckmanEM}}
#'
#' @references
#' Insert any relevant references here.
#'
#' @author Marcos Oliveira
#'
#' @export Influence
#'

Influence <- function(object,type,colx = NULL,k=3.5){

  # types = c("case-weight","scale","response",
  #           "exploratory",
  #           "exploratory2","exploratory3")

  types = c("case-weight","scale","response","exploratory")



  if(is(object,"HeckmanEM")){
    if(object$family == "CN" || object$family == "cn"){
      stop("influence analysis is not available for the contaminated normal model")
    }

    if(!is.numeric(k) | k <= 0){
      stop("Constant k must be a positive real number.")
    }

    if(type %in% types){
      Perturb = which(types == type)
    }else{
      #stop("Perturbation types may only take values `case-weight`,`scale`,`response`,`exploratory`, `exploratory2`, `exploratory3`.")
      stop("Perturbation types may only take values `case-weight`,`scale`,`response`, and `exploratory`.")
    }

    if(Perturb>=4){

      if(is.null(colx)){
        stop("The parameter colx should be an integer that represents the position
         of the column in matrix x that will undergo perturbation.")
      }

      if(!is.numeric(colx) | colx %% 1 != 0 | colx <= 0 | colx > ncol(object$x)){
        stop("The parameter colx should be an integer that represents the position
         of the column in matrix x that will undergo perturbation.")
      }

    }

    If(y = object$y,x = object$x,w = object$w,cc = object$cc,
       beta=object$beta, gama=object$gamma,rho=object$rho,
       sigma=object$sigma,nu = object$nu,
       Perturb = Perturb, j = colx, k=k, family = object$family)

  }else{
    stop("The object must be of the `HeckmanEM` class.")
  }
}


If <- function(y, x, w, cc, beta=beta, gama=gama,rho=rho,
               sigma=sigma, nu, Perturb, j, k=3.5, family = "Normal"){

  # start.time <- Sys.time()

  HQ   <- HessianaQ(y, x, w, cc, beta=beta, gama=gama, rho=rho, sigma=sigma, nu, type = family)

  n        <- nrow(x)
  y        <- matrix(y, n, 1)
  p        <- ncol(x)
  q        <- ncol(w)

  sigma2   <- sigma^2
  rhoa     <- sigma*rho
  Sigma    <- matrix(c(sigma2, rhoa, rhoa, 1), ncol = 2)
  betaC    <- c(beta,gama)

  SigI        <- solve(Sigma) # Inversa de Sigma

  d.Ssigma2     <- matrix(c(1, rho/(2*sigma), rho/(2*sigma), 0), 2, 2)
  d.Srho        <- matrix(c(0, sigma, sigma, 0), 2, 2)

  d2.Ssigma2    <- matrix(c(0, -rho/(4*sigma^3), -rho/(4*sigma^3), 0), 2, 2)
  d2.Ssigma2rho <- matrix(c(0, 1/(2*sigma), 1/(2*sigma), 0), 2, 2)
  d2.Srho       <- matrix(0, 2, 2)

  I2       <- matrix(1, nrow=2)
  delta1   <- c()
  delta2   <- c()
  delta3   <- c()
  mdelta   <- matrix(0, (p+q+2), n)

  if(family=="Normal"){
    mu1     <- x%*%beta
    mu2     <- w%*%gama

    suma1   <- matrix(0, p+q, 1) # Delta_betacw0
    suma2   <- 0                 # Delta_sigma2w0
    suma3   <- 0                 # Delta_rhow0

    if(Perturb==1){ # Case-weight Perturbation
      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          MomNT     <- meanvarTMD(0, Inf, mu12.1, sigma12.1, dist = "normal")

          uy[2]     <- MomNT$mean
          uyy[2,2]  <- MomNT$varcov
          uyy       <- uyy + uy%*%t(uy)
        }
        else{
          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma, dist = "normal")
          uy        <- MomNT1$mean
          uyy       <- MomNT1$EYY
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        aux5        <- aux41%*%betaC                                      # Xic * Bc
        aux2        <- uyy - uy%*%t(aux5) - aux5%*%t(uy) + aux5%*%t(aux5) # Gamma_i

        delta1       <- t(aux41)%*%SigI%*%uy - t(aux41)%*%SigI%*%aux5
        delta2       <- - 0.5*sum(diag(SigI%*%d.Ssigma2)) + 0.5*sum(diag(aux2%*%SigI%*%d.Ssigma2%*%SigI))
        delta3       <- - 0.5*sum(diag(SigI%*%d.Srho))    + 0.5*sum(diag(aux2%*%SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==2){ # Scale Matrix Perturbation
      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          MomNT     <- meanvarTMD(0, Inf, mu12.1, sigma12.1, dist = "normal")

          uy[2]     <- MomNT$mean
          uyy[2,2]  <- MomNT$varcov
          uyy       <- uyy + uy%*%t(uy)
        }
        else{
          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma, dist = "normal")
          uy        <- MomNT1$mean
          uyy       <- MomNT1$EYY
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        aux5        <- aux41%*%betaC                                      # Xic * Bc
        aux2        <- uyy - uy%*%t(aux5) - aux5%*%t(uy) + aux5%*%t(aux5) # Gamma_i

        delta1       <- t(aux41)%*%SigI%*%uy - t(aux41)%*%SigI%*%aux5
        delta2       <- 0.5*sum(diag(aux2%*%SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*sum(diag(aux2%*%SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==3){ # Response Perturbation
      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          MomNT     <- meanvarTMD(0, Inf, mu12.1, sigma12.1, dist = "normal")

          uy[2]     <- MomNT$mean
          uyy[2,2]  <- MomNT$varcov
          uyy       <- uyy + uy%*%t(uy)
        }
        else{
          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma, dist = "normal")
          uy        <- MomNT1$mean
          uyy       <- MomNT1$EYY
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        aux5        <- aux41%*%betaC                                      # Xic * Bc

        a           <- -2*t(uy)%*%I2 + 2*t(aux5)%*%I2

        delta1       <- - t(aux41)%*%SigI%*%I2
        delta2       <- 0.5*a*sum(diag(SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*a*sum(diag(SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==4){ # Explanatory variables Perturbation: Case 1 - Primary Regression
      Iu       <- matrix(0, ncol = p)
      Iu[j]    <- 1
      C        <- matrix(0, nrow = 2, ncol = p+q)
      C[1,1:p] <- Iu

      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          MomNT     <- meanvarTMD(0, Inf, mu12.1, sigma12.1, dist = "normal")

          uy[2]     <- MomNT$mean
          uyy[2,2]  <- MomNT$varcov
          uyy       <- uyy + uy%*%t(uy)
        }
        else{
          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma, dist = "normal")
          uy        <- MomNT1$mean
          uyy       <- MomNT1$EYY
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        #aux41w      <- rbind(c(x[i,]+Iu, 0*w[i,]), c(x[i,]*0, w[i,]))     # Xicw
        aux5        <- aux41%*%betaC                                      # Xic * Bc
        #aux5w       <- aux41w%*%betaC                                     # Xicw * Bc

        aux2        <- -t(uy)%*%C%*%betaC - t(C%*%betaC)%*%uy + t(C%*%betaC)%*%aux5 + t(aux5)%*%C%*%betaC # Gamma_i atualizado

        delta1       <- t(C)%*%SigI%*%uy - t(C)%*%SigI%*%aux5 - t(aux41)%*%SigI%*%C%*%betaC
        delta2       <- 0.5*aux2*sum(diag(SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*aux2*sum(diag(SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==5){ # Explanatory variables Perturbation: Case 2 - Selection equation
      Iu               <- matrix(0, ncol = q)
      Iu[j]            <- 1
      D                <- matrix(0, nrow = 2, ncol = p+q)
      D[2,(p+1):(p+q)] <- Iu

      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          MomNT     <- meanvarTMD(0, Inf, mu12.1, sigma12.1, dist = "normal")

          uy[2]     <- MomNT$mean
          uyy[2,2]  <- MomNT$varcov
          uyy       <- uyy + uy%*%t(uy)
        }
        else{
          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma, dist = "normal")
          uy        <- MomNT1$mean
          uyy       <- MomNT1$EYY
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        #aux41w      <- rbind(c(x[i,]+Iu, 0*w[i,]), c(x[i,]*0, w[i,]))     # Xicw
        aux5        <- aux41%*%betaC                                      # Xic * Bc
        #aux5w       <- aux41w%*%betaC                                     # Xicw * Bc

        aux2        <- -t(uy)%*%D%*%betaC - t(D%*%betaC)%*%uy + t(D%*%betaC)%*%aux5 + t(aux5)%*%D%*%betaC # Gamma_i atualizado

        delta1       <- t(D)%*%SigI%*%uy - t(D)%*%SigI%*%aux5 - t(aux41)%*%SigI%*%D%*%betaC
        delta2       <- 0.5*aux2*sum(diag(SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*aux2*sum(diag(SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==6){ # Explanatory variables Perturbation: Case 3 - Both primary and selection equation
      Iu               <- matrix(0, ncol = p)
      Iv               <- matrix(0, ncol = q)
      Iu[j]            <- 1
      Iv[j]            <- 1
      E                <- matrix(0, nrow = 2, ncol = p+q)
      E[1,1:p]         <- Iu
      E[2,(p+1):(p+q)] <- Iv

      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          MomNT     <- meanvarTMD(0, Inf, mu12.1, sigma12.1, dist = "normal")

          uy[2]     <- MomNT$mean
          uyy[2,2]  <- MomNT$varcov
          uyy       <- uyy + uy%*%t(uy)
        }
        else{
          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma, dist = "normal")
          uy        <- MomNT1$mean
          uyy       <- MomNT1$EYY
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        #aux41w      <- rbind(c(x[i,]+Iu, 0*w[i,]), c(x[i,]*0, w[i,]))     # Xicw
        aux5        <- aux41%*%betaC                                      # Xic * Bc
        #aux5w       <- aux41w%*%betaC                                     # Xicw * Bc

        aux2        <- -t(uy)%*%E%*%betaC - t(E%*%betaC)%*%uy + t(E%*%betaC)%*%aux5 + t(aux5)%*%E%*%betaC # Gamma_i atualizado

        delta1       <- t(E)%*%SigI%*%uy - t(E)%*%SigI%*%aux5 - t(aux41)%*%SigI%*%E%*%betaC
        delta2       <- 0.5*aux2*sum(diag(SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*aux2*sum(diag(SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }
  }

  if(family=="T"){
    mu1     <- x%*%beta
    mu2     <- w%*%gama

    suma1   <- matrix(0, p+q, 1) # Delta_betacw0
    suma2   <- 0                 # Delta_sigma2w0
    suma3   <- 0                 # Delta_rhow0

    if(Perturb==1){ # Case-weight Perturbation
      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          PsiA      <- Sigma*nu/(nu+2)
          nu1       <- (nu + 1)                            # X
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          ScA       <- nu/(nu + 2)*sigma12.1

          Qy1       <- (y[i] - mu1[i])^2/sigma2            # X
          Qy2       <- (y[i] - mu1[i])^2/PsiA[1, 1]        # X
          auxcte    <- as.numeric((nu + Qy1)/(nu + 1))
          auxcte1   <- as.numeric((nu + 2 + Qy2)/(nu + 3))

          Sc22      <- auxcte*sigma12.1

          muUi      <- mu12.1
          SigmaUi   <- Sc22

          SigmaUiA  <- auxcte1*ScA
          auxupper  <- -muUi

          auxU1     <- 1 - pt(auxupper/sqrt(SigmaUiA), nu1 + 2)
          auxU2     <- 1 - pt(auxupper/sqrt(SigmaUi), nu1)

          MoMT      <- meanvarTMD(0, Inf, muUi, SigmaUiA, dist = "t", nu = nu1 + 2)

          vary      <- matrix(0, 2, 2)
          vary[2,2] <- meanvarTMD(0, Inf, muUi, SigmaUi , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)/auxcte
          U1        <- (U0)*(MoMT$mean)
          U2        <- (U0)*(MoMT$EYY)

          Auxtuy    <- (matrix(y[i], 2, 1))

          uy        <- Auxtuy*U0
          uy[2]     <- U1

          uyy       <- (Auxtuy%*%t(Auxtuy))

          AAx       <- uyy[1, 1]*U0
          ABx       <- Auxtuy[1]%*%t(U1)
          BBx       <- U2

          uyy[1,1]  <- AAx
          uyy[1,2]  <- ABx
          uyy[2,1]  <- ABx
          uyy[2,2]  <- BBx
        }
        else{
          SigmaUi   <- Sigma
          SigmaUiA  <- Sigma*nu/(nu+2)
          auxupper  <- -mu2[i]

          auxU1     <- pt(auxupper/sqrt(SigmaUiA[2,2]), nu+2)
          auxU2     <- pt(auxupper/sqrt(SigmaUi[2,2]) , nu)

          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), SigmaUiA, dist = "t", nu = nu+2)
          vary      <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma   , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)
          U1        <- auxU1/auxU2*MomNT1$mean
          U2        <- auxU1/auxU2*MomNT1$EYY

          uy        <- U1
          uyy       <- U2
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))           # Xic
        aux5        <- aux41%*%betaC                                         # Xic * Bc
        aux2        <- uyy - uy%*%t(aux5) - aux5%*%t(uy) + U0*aux5%*%t(aux5) # Gamma_i atualizado

        delta1       <- t(aux41)%*%SigI%*%uy - U0*t(aux41)%*%SigI%*%aux5
        delta2       <- - 0.5*sum(diag(SigI%*%d.Ssigma2)) + 0.5*sum(diag(aux2%*%SigI%*%d.Ssigma2%*%SigI))
        delta3       <- - 0.5*sum(diag(SigI%*%d.Srho))    + 0.5*sum(diag(aux2%*%SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==2){ # Scale Matrix Perturbation
      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          PsiA      <- Sigma*nu/(nu+2)
          nu1       <- (nu + 1)                            # X
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          ScA       <- nu/(nu + 2)*sigma12.1

          Qy1       <- (y[i] - mu1[i])^2/sigma2            # X
          Qy2       <- (y[i] - mu1[i])^2/PsiA[1, 1]        # X
          auxcte    <- as.numeric((nu + Qy1)/(nu + 1))
          auxcte1   <- as.numeric((nu + 2 + Qy2)/(nu + 3))

          Sc22      <- auxcte*sigma12.1

          muUi      <- mu12.1
          SigmaUi   <- Sc22

          SigmaUiA  <- auxcte1*ScA
          auxupper  <- -muUi

          auxU1     <- 1 - pt(auxupper/sqrt(SigmaUiA), nu1 + 2)
          auxU2     <- 1 - pt(auxupper/sqrt(SigmaUi), nu1)

          MoMT      <- meanvarTMD(0, Inf, muUi, SigmaUiA, dist = "t", nu = nu1 + 2)

          vary      <- matrix(0, 2, 2)
          vary[2,2] <- meanvarTMD(0, Inf, muUi, SigmaUi , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)/auxcte
          U1        <- (U0)*(MoMT$mean)
          U2        <- (U0)*(MoMT$EYY)

          Auxtuy    <- (matrix(y[i], 2, 1))

          uy        <- Auxtuy*U0
          uy[2]     <- U1

          uyy       <- (Auxtuy%*%t(Auxtuy))

          AAx       <- uyy[1, 1]*U0
          ABx       <- Auxtuy[1]%*%t(U1)
          BBx       <- U2

          uyy[1,1]  <- AAx
          uyy[1,2]  <- ABx
          uyy[2,1]  <- ABx
          uyy[2,2]  <- BBx
        }
        else{
          SigmaUi   <- Sigma
          SigmaUiA  <- Sigma*nu/(nu+2)
          auxupper  <- -mu2[i]

          auxU1     <- pt(auxupper/sqrt(SigmaUiA[2,2]), nu+2)
          auxU2     <- pt(auxupper/sqrt(SigmaUi[2,2]) , nu)

          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), SigmaUiA, dist = "t", nu = nu+2)
          vary      <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma   , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)
          U1        <- auxU1/auxU2*MomNT1$mean
          U2        <- auxU1/auxU2*MomNT1$EYY

          uy        <- U1
          uyy       <- U2
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))           # Xic
        aux5        <- aux41%*%betaC                                         # Xic * Bc
        aux2        <- uyy - uy%*%t(aux5) - aux5%*%t(uy) + U0*aux5%*%t(aux5) # Gamma_i atualizado

        delta1       <- t(aux41)%*%SigI%*%uy - U0*t(aux41)%*%SigI%*%aux5
        delta2       <- 0.5*sum(diag(aux2%*%SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*sum(diag(aux2%*%SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==3){ # Response Perturbation
      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          PsiA      <- Sigma*nu/(nu+2)
          nu1       <- (nu + 1)                            # X
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          ScA       <- nu/(nu + 2)*sigma12.1

          Qy1       <- (y[i] - mu1[i])^2/sigma2            # X
          Qy2       <- (y[i] - mu1[i])^2/PsiA[1, 1]        # X
          auxcte    <- as.numeric((nu + Qy1)/(nu + 1))
          auxcte1   <- as.numeric((nu + 2 + Qy2)/(nu + 3))

          Sc22      <- auxcte*sigma12.1

          muUi      <- mu12.1
          SigmaUi   <- Sc22

          SigmaUiA  <- auxcte1*ScA
          auxupper  <- -muUi

          auxU1     <- 1 - pt(auxupper/sqrt(SigmaUiA), nu1 + 2)
          auxU2     <- 1 - pt(auxupper/sqrt(SigmaUi), nu1)

          MoMT      <- meanvarTMD(0, Inf, muUi, SigmaUiA, dist = "t", nu = nu1 + 2)

          vary      <- matrix(0, 2, 2)
          vary[2,2] <- meanvarTMD(0, Inf, muUi, SigmaUi , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)/auxcte
          U1        <- (U0)*(MoMT$mean)
          U2        <- (U0)*(MoMT$EYY)

          Auxtuy    <- (matrix(y[i], 2, 1))

          uy        <- Auxtuy*U0
          uy[2]     <- U1

          uyy       <- (Auxtuy%*%t(Auxtuy))

          AAx       <- uyy[1, 1]*U0
          ABx       <- Auxtuy[1]%*%t(U1)
          BBx       <- U2

          uyy[1,1]  <- AAx
          uyy[1,2]  <- ABx
          uyy[2,1]  <- ABx
          uyy[2,2]  <- BBx
        }
        else{
          SigmaUi   <- Sigma
          SigmaUiA  <- Sigma*nu/(nu+2)
          auxupper  <- -mu2[i]

          auxU1     <- pt(auxupper/sqrt(SigmaUiA[2,2]), nu+2)
          auxU2     <- pt(auxupper/sqrt(SigmaUi[2,2]) , nu)

          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), SigmaUiA, dist = "t", nu = nu+2)
          vary      <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma   , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)
          U1        <- auxU1/auxU2*MomNT1$mean
          U2        <- auxU1/auxU2*MomNT1$EYY

          uy        <- U1
          uyy       <- U2
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))           # Xic
        aux5        <- aux41%*%betaC                                         # Xic * Bc

        a           <- -2*t(uy)%*%I2 + 2*U0*t(aux5)%*%I2

        delta1       <- - U0*t(aux41)%*%SigI%*%I2
        delta2       <- 0.5*a*sum(diag(SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*a*sum(diag(SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==4){ # Explanatory variables Perturbation: Case 1 - Primary Regression
      Iu       <- matrix(0, ncol = p)
      Iu[j]    <- 1
      C        <- matrix(0, nrow = 2, ncol = p+q)
      C[1,1:p] <- Iu

      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          PsiA      <- Sigma*nu/(nu+2)
          nu1       <- (nu + 1)                            # X
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          ScA       <- nu/(nu + 2)*sigma12.1

          Qy1       <- (y[i] - mu1[i])^2/sigma2            # X
          Qy2       <- (y[i] - mu1[i])^2/PsiA[1, 1]        # X
          auxcte    <- as.numeric((nu + Qy1)/(nu + 1))
          auxcte1   <- as.numeric((nu + 2 + Qy2)/(nu + 3))

          Sc22      <- auxcte*sigma12.1

          muUi      <- mu12.1
          SigmaUi   <- Sc22

          SigmaUiA  <- auxcte1*ScA
          auxupper  <- -muUi

          auxU1     <- 1 - pt(auxupper/sqrt(SigmaUiA), nu1 + 2)
          auxU2     <- 1 - pt(auxupper/sqrt(SigmaUi), nu1)

          MoMT      <- meanvarTMD(0, Inf, muUi, SigmaUiA, dist = "t", nu = nu1 + 2)

          vary      <- matrix(0, 2, 2)
          vary[2,2] <- meanvarTMD(0, Inf, muUi, SigmaUi , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)/auxcte
          U1        <- (U0)*(MoMT$mean)
          U2        <- (U0)*(MoMT$EYY)

          Auxtuy    <- (matrix(y[i], 2, 1))

          uy        <- Auxtuy*U0
          uy[2]     <- U1

          uyy       <- (Auxtuy%*%t(Auxtuy))

          AAx       <- uyy[1, 1]*U0
          ABx       <- Auxtuy[1]%*%t(U1)
          BBx       <- U2

          uyy[1,1]  <- AAx
          uyy[1,2]  <- ABx
          uyy[2,1]  <- ABx
          uyy[2,2]  <- BBx
        }
        else{
          SigmaUi   <- Sigma
          SigmaUiA  <- Sigma*nu/(nu+2)
          auxupper  <- -mu2[i]

          auxU1     <- pt(auxupper/sqrt(SigmaUiA[2,2]), nu+2)
          auxU2     <- pt(auxupper/sqrt(SigmaUi[2,2]) , nu)

          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), SigmaUiA, dist = "t", nu = nu+2)
          vary      <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma   , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)
          U1        <- auxU1/auxU2*MomNT1$mean
          U2        <- auxU1/auxU2*MomNT1$EYY

          uy        <- U1
          uyy       <- U2
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        #aux41w      <- rbind(c(x[i,]+Iu, 0*w[i,]), c(x[i,]*0, w[i,]))     # Xicw
        aux5        <- aux41%*%betaC                                      # Xic * Bc
        #aux5w       <- aux41w%*%betaC                                     # Xicw * Bc

        aux2        <- -t(uy)%*%C%*%betaC - t(C%*%betaC)%*%uy + U0*t(C%*%betaC)%*%aux5 + U0*t(aux5)%*%C%*%betaC # Gamma_i atualizado

        delta1       <- t(C)%*%SigI%*%uy - U0*t(C)%*%SigI%*%aux5 - U0*t(aux41)%*%SigI%*%C%*%betaC
        delta2       <- 0.5*aux2*sum(diag(SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*aux2*sum(diag(SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==5){ # Explanatory variables Perturbation: Case 2 - Selection equation
      Iu               <- matrix(0, ncol = q)
      Iu[j]            <- 1
      D                <- matrix(0, nrow = 2, ncol = p+q)
      D[2,(p+1):(p+q)] <- Iu

      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          PsiA      <- Sigma*nu/(nu+2)
          nu1       <- (nu + 1)                            # X
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          ScA       <- nu/(nu + 2)*sigma12.1

          Qy1       <- (y[i] - mu1[i])^2/sigma2            # X
          Qy2       <- (y[i] - mu1[i])^2/PsiA[1, 1]        # X
          auxcte    <- as.numeric((nu + Qy1)/(nu + 1))
          auxcte1   <- as.numeric((nu + 2 + Qy2)/(nu + 3))

          Sc22      <- auxcte*sigma12.1

          muUi      <- mu12.1
          SigmaUi   <- Sc22

          SigmaUiA  <- auxcte1*ScA
          auxupper  <- -muUi

          auxU1     <- 1 - pt(auxupper/sqrt(SigmaUiA), nu1 + 2)
          auxU2     <- 1 - pt(auxupper/sqrt(SigmaUi), nu1)

          MoMT      <- meanvarTMD(0, Inf, muUi, SigmaUiA, dist = "t", nu = nu1 + 2)

          vary      <- matrix(0, 2, 2)
          vary[2,2] <- meanvarTMD(0, Inf, muUi, SigmaUi , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)/auxcte
          U1        <- (U0)*(MoMT$mean)
          U2        <- (U0)*(MoMT$EYY)

          Auxtuy    <- (matrix(y[i], 2, 1))

          uy        <- Auxtuy*U0
          uy[2]     <- U1

          uyy       <- (Auxtuy%*%t(Auxtuy))

          AAx       <- uyy[1, 1]*U0
          ABx       <- Auxtuy[1]%*%t(U1)
          BBx       <- U2

          uyy[1,1]  <- AAx
          uyy[1,2]  <- ABx
          uyy[2,1]  <- ABx
          uyy[2,2]  <- BBx
        }
        else{
          SigmaUi   <- Sigma
          SigmaUiA  <- Sigma*nu/(nu+2)
          auxupper  <- -mu2[i]

          auxU1     <- pt(auxupper/sqrt(SigmaUiA[2,2]), nu+2)
          auxU2     <- pt(auxupper/sqrt(SigmaUi[2,2]) , nu)

          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), SigmaUiA, dist = "t", nu = nu+2)
          vary      <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma   , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)
          U1        <- auxU1/auxU2*MomNT1$mean
          U2        <- auxU1/auxU2*MomNT1$EYY

          uy        <- U1
          uyy       <- U2
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        #aux41w      <- rbind(c(x[i,]+Iu, 0*w[i,]), c(x[i,]*0, w[i,]))     # Xicw
        aux5        <- aux41%*%betaC                                      # Xic * Bc
        #aux5w       <- aux41w%*%betaC                                     # Xicw * Bc

        aux2        <- -t(uy)%*%D%*%betaC - t(D%*%betaC)%*%uy + U0*t(D%*%betaC)%*%aux5 + U0*t(aux5)%*%D%*%betaC # Gamma_i atualizado

        delta1       <- t(D)%*%SigI%*%uy - U0*t(D)%*%SigI%*%aux5 - U0*t(aux41)%*%SigI%*%D%*%betaC
        delta2       <- 0.5*aux2*sum(diag(SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*aux2*sum(diag(SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }

    if(Perturb==6){ # Explanatory variables Perturbation: Case 6 - Both primary and selection equation
      Iu               <- matrix(0, ncol = p)
      Iv               <- matrix(0, ncol = q)
      Iu[j]            <- 1
      Iv[j]            <- 1
      E                <- matrix(0, nrow = 2, ncol = p+q)
      E[1,1:p]         <- Iu
      E[2,(p+1):(p+q)] <- Iv

      for (i in 1:n){
        uy  <- matrix(y[i],2,1)
        uyy <- matrix(0,2,2)

        if(cc[i]==1){
          PsiA      <- Sigma*nu/(nu+2)
          nu1       <- (nu + 1)                            # X
          mu12.1    <- mu2[i] + rho/sigma*(y[i] - mu1[i])
          sigma12.1 <- 1 - rho^2
          ScA       <- nu/(nu + 2)*sigma12.1

          Qy1       <- (y[i] - mu1[i])^2/sigma2            # X
          Qy2       <- (y[i] - mu1[i])^2/PsiA[1, 1]        # X
          auxcte    <- as.numeric((nu + Qy1)/(nu + 1))
          auxcte1   <- as.numeric((nu + 2 + Qy2)/(nu + 3))

          Sc22      <- auxcte*sigma12.1

          muUi      <- mu12.1
          SigmaUi   <- Sc22

          SigmaUiA  <- auxcte1*ScA
          auxupper  <- -muUi

          auxU1     <- 1 - pt(auxupper/sqrt(SigmaUiA), nu1 + 2)
          auxU2     <- 1 - pt(auxupper/sqrt(SigmaUi), nu1)

          MoMT      <- meanvarTMD(0, Inf, muUi, SigmaUiA, dist = "t", nu = nu1 + 2)

          vary      <- matrix(0, 2, 2)
          vary[2,2] <- meanvarTMD(0, Inf, muUi, SigmaUi , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)/auxcte
          U1        <- (U0)*(MoMT$mean)
          U2        <- (U0)*(MoMT$EYY)

          Auxtuy    <- (matrix(y[i], 2, 1))

          uy        <- Auxtuy*U0
          uy[2]     <- U1

          uyy       <- (Auxtuy%*%t(Auxtuy))

          AAx       <- uyy[1, 1]*U0
          ABx       <- Auxtuy[1]%*%t(U1)
          BBx       <- U2

          uyy[1,1]  <- AAx
          uyy[1,2]  <- ABx
          uyy[2,1]  <- ABx
          uyy[2,2]  <- BBx
        }
        else{
          SigmaUi   <- Sigma
          SigmaUiA  <- Sigma*nu/(nu+2)
          auxupper  <- -mu2[i]

          auxU1     <- pt(auxupper/sqrt(SigmaUiA[2,2]), nu+2)
          auxU2     <- pt(auxupper/sqrt(SigmaUi[2,2]) , nu)

          MomNT1    <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), SigmaUiA, dist = "t", nu = nu+2)
          vary      <- meanvarTMD(c(-Inf, -Inf), c(Inf, 0), c(mu1[i], mu2[i]), Sigma   , dist = "t", nu = nu)$varcov

          U0        <- as.numeric(auxU1/auxU2)
          U1        <- auxU1/auxU2*MomNT1$mean
          U2        <- auxU1/auxU2*MomNT1$EYY

          uy        <- U1
          uyy       <- U2
        }

        aux41       <- rbind(c(x[i,], 0*w[i,]), c(x[i,]*0, w[i,]))        # Xic
        #aux41w      <- rbind(c(x[i,]+Iu, 0*w[i,]), c(x[i,]*0, w[i,]))     # Xicw
        aux5        <- aux41%*%betaC                                      # Xic * Bc
        #aux5w       <- aux41w%*%betaC                                     # Xicw * Bc

        aux2        <- -t(uy)%*%E%*%betaC - t(E%*%betaC)%*%uy + U0*t(E%*%betaC)%*%aux5 + U0*t(aux5)%*%E%*%betaC # Gamma_i atualizado

        delta1       <- t(E)%*%SigI%*%uy - U0*t(E)%*%SigI%*%aux5 - U0*t(aux41)%*%SigI%*%E%*%betaC
        delta2       <- 0.5*aux2*sum(diag(SigI%*%d.Ssigma2%*%SigI))
        delta3       <- 0.5*aux2*sum(diag(SigI%*%d.Srho%*%SigI))

        mdelta[1:(p+q),i] <- delta1
        mdelta[p+q+1,i]   <- delta2
        mdelta[p+q+2,i]   <- delta3
      }
    }
  }

  If        <- t(mdelta)%*%solve(HQ$MatrizQ)%*%mdelta
  medida    <- abs(diag(If)/sum(diag(If))) # Verificar se está correto utilizar "abs", para pegar somente pontos positivos
  benchmark <- mean(medida) + k*sd(medida)

  obs <- c(rep(0,n))

  for(i in 1:n){
    if(medida[i] > benchmark){ obs[i] <- i }
  }

  influent <- obs[obs!=0]

  if(length(influent) == 0){influent <- NULL}

  obj.out <- list(M0 = medida, benchmark = benchmark,
                  influent = influent, type = c("case-weight","scale","response","exploratory")[Perturb])

  class(obj.out) = "HeckmanEM.influence"

  return(obj.out)


}


# -------------------------------------------------------------------------
globalVariables(c("GD","M0"))
#' @importFrom methods is
#' @importFrom stats sd
#' @importFrom ggplot2 ggplot geom_segment theme_bw geom_hline aes labs
#' @export
plot.HeckmanEM.deletion = function(x, ...) {
  GDplot = data.frame(x=seq(1, length(c(x$GD))), GD=x$GD)
  m = nrow(GDplot)
  bench = x$benchmark

  ggplot(GDplot) + geom_segment(aes(x=x, xend=x, y=0, yend=GD)) + theme_bw() +
    labs(x="Index", y="GD") + geom_hline(yintercept=bench, color="red", linetype="dotted")
}


# -------------------------------------------------------------------------

#' @export
plot.HeckmanEM.influence = function(x, ...) {
  M0plot = data.frame(x=seq(1, length(c(x$M0))), M0=x$M0)
  m = nrow(M0plot)
  bench = x$benchmark
  type = x$type

  ggplot(M0plot, aes(x=x, y=M0)) + geom_segment(aes(x=x, xend=x, y=0, yend=M0)) +
    geom_hline(yintercept=x$benchmark, color="red", linetype="dotted") +
    labs(x="Index", y=bquote(M[0]),
         subtitle=paste("Perturbation type:",type)) + theme_bw()
}
