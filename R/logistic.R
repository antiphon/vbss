#' Fit VB Spike and Slab logistic regression
#'
#' The workhorse for logistic \link{vbglmss}.
#' 
#' 
#' @param y binary response.
#' @param X the SS covariates.
#' @param Z the non-SS (regular Gaussian) covariates.
#' @param verbose print additional runtime messages
#' @param eps the approximation loop convergence threshold
#' @param ... further parameters passed to sub-functions.
#' This function makes the decision of which sub-functions need to be called
#' based on the presense/absence of X and Z and parameter values.
#' @import stringr Matrix
#' @export
vbglmss.fit.logistic <- function(y, X=NULL, Z=NULL, 
                                 prior=list(pi=0.5, tau=list(a=10,b=10)), verbose=TRUE,
                                 ...){
  ##
  ## check prior
  init <- list(mw=T, mb=T, tau=T, pi=T)
  ##
  ## verbosity
  cat2 <- function(x)NULL
  if(verbose) cat2<-cat
  cat2("Logistic family.\n")
  ##
  ## we use -1/+1 formulation
  y[y<1] <- -1
  ##
  ## check coefficient models
  ##
  ## The gaussian coef's:
  gaussian <- !is.null(Z)
  if(gaussian){
    ## we use transposes
    Z <- t(Z)
    P  <- nrow(Z)
    mb <- Matrix(rep(0, P))
    C <- Z%*%y
    mb0 <- mb
    Sbi0 <- Sb0  <- Diagonal(n=P)
    cat2(P," gaussian coefficients.\n")
  }
  ##
  ## Then the SS coef's:
  nongaussian <- !is.null(X)
  if(nongaussian){
    K <- ncol(X)
    prior$SS <- list()
    prior$SS$tau <- prior$tau
    if(init$mw) prior$SS$mw <- rep(0, K)
    if(init$pi) prior$SS$pi <- rep(prior$pi, K)
    if(init$tau) prior$SS$tau <- rep(1, K)
    cat2(K," spike-n-slab coefficients.\n")
  }
  
  ## call
  if(!gaussian & nongaussian) f<-vbglmss.fit.logistic.X
  else if(gaussian & !nongaussian) f<-vbglmss.fit.logistic.Z
  else if(gaussian & nongaussian) f<-vbglmss.fit.logistic.XZ
  
  f(X=X, y=y, Z=Z, prior=prior, cat2=cat2, ...)
### End of VB logistic regression parent function.
}


#' Fit logistic regression with only SS coefficients.
#' 
#' @param y binary response.
#' @param X the SS covariates.
#' @param prior a list of priors
#' @param cat2 printer
#' @param eps the approximation loop convergence threshold
#' @export
vbglmss.fit.logistic.X <- function(y, X, prior, cat2, eps=1e-2, ...){
  N <- length(y)
  X <- t(Matrix(X))
  K <- nrow(X)
  D <- X%*%y   
  m0 <- prior$SS$m
  m <- m0
  tau <- prior$SS$tau
  pi0 <- prior$SS$pi
  g <- pi0
  xi <- rep(5, N)
  Et <- 1
  logp<- -Inf
  logp_hist <- log_const <- 0
  gprior <- log(pi0/(1-pi0))
  ##################################
  ## The optimization loop:
  loop <- TRUE
  while(loop){
    ## collect for convergence diagnostics:
    gold <- g
    xiold <- xi
    ## update the logistic link approximation parameters
    la <- lambda(as.numeric(xi))
    L <- Diagonal(x=la)
    B <- X%*%L%*%t(X)
    ##
    tau <- Et + 2*diag(B)
    v <- (0.5*D + m0*Et)/tau
    for(k in 1:K){
      e      <- sum(g[-k]*m[-k]*B[k,-k])
      m[k]   <- v[k] - 2*e/tau[k]
    }
    u   <- 0.5*(m^2*tau - Et*m0^2) + gprior
    g   <- logit(u)
    ## hyper tau    
    w <- g*m
    w2 <- g*( 1/tau+m^2 )
    Et <- ( a<-(prior$tau$a + K/2)  )/( b<- (prior$tau$b + 0.5*sum(w2 - 2*w*m0 + m0^2)) )
    ## variational parameters
    xi2 <- diag( t(X)%*%(Diagonal(x=1/tau)+m%*%t(m))%*%X )
    xi <- sqrt(as.numeric(xi2))
    ##
    ## marginal likelihood
    logp <-  
    ##
    d <- max(abs(xi-xiold), abs(g-gold))
    loop<- d > eps
    cat2("\r",d, "         ")
  }
  logp <- logp_const + logp_hist
  list(SS=list(m=m, tau=tau, g=g, xi=xi, hypertau=c(a=a,b=b)), logp=logp)
}



#' Fit logistic with only gaussian coef covariates 
#' 
#' @param y binary response.
#' @param Z the covariates.
#' @param prior a list of priors
#' @param cat2 printer
#' @param eps the approximation loop convergence threshold
#' @export

vbglmss.fit.logistic.Z <- function(y, Z, prior, cat2, eps, ...){
  stop("logistic.Z not implemented.")
  N <- length(y)
  X <- t(X)
  K <- nrow(X)
  D <- X%*%y   
  m0 <- prior$mw
  m <- m0
  tau <- prior$tau
  pi <- prior$pi
  xi <- 
    ##################################
  ## The optimization loop:
  loop <- TRUE
  while(loop){
    ## collect for convergence diagnostics:
    piold <- pi
    xiold <- xi
    ## update the logistic link approximation parameters
    la <- lambda(as.numeric(xi))
    ##
  }
}

vbglmss.fit.logistic.XZ <- function(y, X, Z, prior, cat2, eps, ...){
  stop("logistic.XZ not implemented.")
  N <- length(y)
  X <- t(X)
  K <- nrow(X)
  D <- X%*%y   
  m0 <- prior$mw
  m <- m0
  tau <- prior$tau
  pi <- prior$pi
  xi <- 
    ##################################
  ## The optimization loop:
  loop <- TRUE
  while(loop){
    ## collect for convergence diagnostics:
    piold <- pi
    xiold <- xi
    ## update the logistic link approximation parameters
    la <- lambda(as.numeric(xi))
    ##
  }
}
