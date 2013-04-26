#' Variational approximation of Spike and Slab regression
#'
#' The function fits a regression model where observations can be from
#' logistic and/or gaussian likelihood. The coefficients are shared, and can have
#'  Spike and Slab (SS) and/or Gaussian (G) priors. 
#' @param y Binary observations
#' @param Xy SS covariates for y
#' @param Zy G covariates for y
#' @param h Gaussian observations
#' @param Xh SS covariates for h
#' @param Zh G covariates for h
#' @param prior A list of priors. See details.
#' @param verbose print additional runtime messages
#' @param damp Damping coefficient in SS mean update. damp * old + (1-damp)*new.
#' @param eps the approximation loop convergence threshold
#' @param max.iter --
#' @details
#' TODO
#' @return 
#'  A list with posteriors.
#' @import Matrix
#' @export

vbglmss.mixedSS<-function(y, Xy, Zy, h, Xh, Zh, 
                     prior=list(theta=list(),
                                beta=list(),
                                binary=list(),
                                gaussian=list()
                     ),
                     verbose=FALSE,
                     damp=0,
                     eps=1e-2, max.iter=100) {
  ## denote by h the gaussian data
  ## and y the binary data.
  ## X are the SS covariates and
  ## Z are the regular covariates.
  ##
  ## helper functions:
  logit <- function(x)1/(1+exp(-x))
  lambda <- function(l)  (logit(l)-0.5)/(2*l)
  cat2 <- function(...)NULL
  if(verbose)cat2<-cat
  ##
  ## check data.
  if(missing(y)&missing(h)) stop("Missing both y and h.")
  Ny <- Nh <- K <- J <- Jy <- Jh <- xi <- 0
  ## initialize generic.
  ## check binary data and covariates
  if(!missing(y)) {
    Ny <- length(y)
    if(missing(Xy) & missing(Zy)) stop("No covariates for binary observations given.")
    if(!missing(Xy)){
      Xy <- Matrix(Xy)
      K <- Ky  <- ncol(Xy)
      if(nrow(Xy) != Ny) stop("Dimension mismatch: y and Xy")
    }
    if(!missing(Zy)) {
      Zy <- Matrix(Zy)
      if(nrow(Zy)!= Ny) stop("Dimension mismatch: y and Zy")
      J <- Jy <- ncol(Zy)
      J
    }
    ## use -1/+1 for y's
    y[y<1] <- -1
    if(is.null(prior$binary$xi))
      xi <- rep(10, Ny)
    else xi <- rep(prior$binary$xi, Ny)[1:Ny]
  }
  # check gaussian data and covariates
  if(!missing(h)) {
    Nh <- length(h)
    if(missing(Xh) & missing(Zh))stop("No covariates for gaussian observations given.")
    if(!missing(Xh)) {
      Xh <- Matrix(Xh)
      K <- Kh  <- ncol(Xh)
      if(Ny) if(Kh != Ky) stop("Dimension mismatch: Xy and Xh")
      if(nrow(Xh)!= Nh) stop("Dimension mismatch: h and Xh")
      tXhXh <- t(Xh)%*%Xh
    }
    if(!missing(Zh)) {
      Zh <- Matrix(Zh)
      if(nrow(Zh)!= Nh) stop("Covariate matrices Xh and Zh differ in sample size")
      J <- Jh <- ncol(Zh)
      if(Jy) if(Jh != Jy) stop("Dimension mismatch: Zy and Zh")
      if(nrow(Zh)!= Nh) stop("Dimension mismatch: h and Zh")
      tZhZh <- t(Zh)%*%Zh      
    }
    tauh_a0 <- c(prior$gaussian$tau[1], 0.001)[1]
    tauh_b0 <- c(prior$gaussian$tau[2], 0.001)[1]
    tauh <- tauh_a0/tauh_b0
  }
  ## damping for mean of slab
  damp1 <- damp
  damp2 <- 1-damp1
  ## initialize
  if(K){
    m0theta <- mmtheta <- mtheta <- rep(0, K)
    tautheta <- rep(1, K)
    if(is.null(prior$theta$pi)) pi0 <- rep(0.01, K)
    else pi0 <- rep(prior$theta$pi, K)[1:K]
    gtheta <- pi0
    
    tauhyper_a0 <- c(prior$theta$tau[1], 0.001)[1]
    tauhyper_b0 <- c(prior$theta$tau[2], 0.001)[1]
    tauhyper <- tauhyper_a0/tauhyper_b0
    
  } else gtheta <- -1
  if(J){
    m0beta <- mbeta <- rep(0, J)
    if(is.null(m0beta<-prior$beta$m)) m0beta <- rep(0, J)
    if(is.null(S0beta<-prior$beta$S))S0beta <- Diagonal(x=rep(1e4, J))
    if(is.null(Si0beta<-prior$beta$Si))Si0beta <- solve(S0beta)
    mbeta <- m0beta
    if(is.null(prior$beta$diagonal)) prior$beta$diagonal <- FALSE
  }else mbeta <- 0
  
  A1 <- A2 <- B1 <- B2 <- 0
  d1 <- d2 <- d3 <- d4 <- e1 <- e2 <- e3 <- e4 <- e5 <- 0
  ##
  cat2("Data:",Ny, Nh,", coefs:",K, J,"\n")
  ##########
  ## run.
  iter <- 0
  loop <- TRUE
  while(loop){ 
    iter <- iter + 1
    xiold <- xi
    mbetaold <- mbeta
    gold <- gtheta
    if(Ny) LAMBDA <- Diagonal(n=Ny, x=lambda(as.numeric(xi)) )
    ## update SS
    if(K){
      ## update slabs
      if(Ny) A1 <- t(Xy)%*%LAMBDA%*%Xy
      if(Nh) A2 <- tauh*tXhXh
      A3 <- 2*A1 + A2
      if(Ny & J) d1 <- t(Xy)%*%LAMBDA%*%Zy%*%mbeta
      if(Nh & J) d2 <- tauh*t(Xh)%*%Zh%*%mbeta
      if(Ny) d3 <- t(Xy)%*%y
      if(Nh) d4 <- tauh*t(Xh)%*%h
      d5 <- 0.5*d3 + d4 - d2 - 2*d1
      tautheta <- diag(A3) + tauhyper
      diag(A3) <- 0
      ## damping on the mean update
      mtheta <- damp1*mtheta+
        damp2*(d5 - A3%*%mmtheta + tauhyper*m0theta)/tautheta
      ## update spikes
      u <- 0.5*mtheta^2*tautheta - 0.5*m0theta^2*tauhyper + log(pi0/(1-pi0))
      gtheta <- as.numeric(logit(u))
      mmtheta <- gtheta*mtheta
      
      Stheta <- Diagonal(x=as.numeric(gtheta/tautheta))
      ## update hypertau
      tauhyper_a <- tauhyper_a0 + K/2  
      tauhyper_b <- tauhyper_b0 + 0.5*sum(tautheta + sum((mmtheta-m0theta)^2) )
      tauhyper <- tauhyper_a/tauhyper_b
    }
    
    ## update gaussian coefs
    if(J){
      if(Ny) B1 <- t(Zy)%*%LAMBDA%*%Zy
      if(Nh) B2 <- tauh*t(Zh)%*%Zh
      B3 <- 2*B1 + B2 + Si0beta
      if(Nh & K) e1 <- tauh*t(Zh)%*%Xh%*%mmtheta
      if(Ny & K) e2 <- t(Zy)%*%LAMBDA%*%Xy%*%mmtheta      
      if(Ny) e3 <- t(Zy)%*%y      
      if(Nh) e4 <- tauh*t(Zh)%*%h
      e5 <- Si0beta%*%m0beta
      e6 <-0.5*e3 + e4 + e5 - e1 - 2*e2
      if(prior$beta$diagonal) B3 <- Diagonal(n=J, x=diag(B3))
      Sbeta <- solve(B3)
      mbeta <- as.numeric(Sbeta%*%e6)
    }
    ## update gaussian observations' precision
    if(Nh){
      s <- 0
      if(K) s <- s + sum(diag(tXhXh %*%Stheta))+t(mmtheta)%*%tXhXh%*%mmtheta +
        t(Xh%*%mmtheta)%*%h
      if(J) s <- s + sum(diag(tZhZh %*% Sbeta))+t(mbeta) %*%tZhZh %*%mbeta +
        t(Zh%*%mbeta)%*%h
      if(K & J) s <- s + 2*t(mbeta)%*%t(Zh)%*%Xh%*%mmtheta
      tauh_a <- as.numeric(tauh_a0 + Nh/2)
      tauh_b <- as.numeric(tauh_b0 + 0.5*s)
      tauh <- tauh_a/tauh_b
    }
    ## update vb logistic approx parameters
    if(Ny){
      w <- rep(0, Ny)
      S <- Diagonal(n=Ny, x=0)
      if(K) S <- S + Xy%*%Stheta%*%t(Xy)
      if(J) S <- S + Zy%*%Sbeta%*%t(Zy)
      if(K) w <- w + Xy%*%mmtheta
      if(J) w <- w + Zy%*%mbeta
      xi <- as.numeric(sqrt( diag(S+w%*%t(w)) ))
    }
    ## loop ending...
    loop<- (iter < max.iter) &
      (d<-max(c(abs(xi-xiold), abs(gold-gtheta), abs(mbeta-mbetaold)))) > eps
    
    cat2("d:", d, "           \r")
  } ## loop end.
  cat2("\n")
  if(iter >= max.iter) warning("Did not converge.")
  ## compile result.
  res <- list(iter=iter, call=sys.call())
  if(K) res$theta <- list(m=mtheta, S=Diagonal(x=1/tautheta), g=gtheta,
                          tauhyper=c(tauhyper_a, tauhyper_b), mmtheta=mmtheta)
  if(J) res$beta <- list(m=mbeta, S=Sbeta)
  if(Ny) res$binary=list(xi=xi)
  if(Nh) res$gaussian <- list(tauh=c(tauh_a, tauh_b))
  ## 
  res
}
### End of VB function.
