#' Fit VB Spike and Slab linear regression
#'
#' The workhorse for linear (Gaussian family) \link{vbglmss}.
#' 
#' 
#' @param y response.
#' @param X the SS covariates
#' @param Z the non-SS (regular Gaussian) covariates
#' @param method Use 'simultaneous' or 'sequential'. See Details.
#' @param verbose print additional runtime messages
#' @param diagonal use diagonal covariance matrix for regular coefficients?
#' @param eps the approximation loop convergence threshold

vbglmss.fit.linear <- function(y, X, Z, method="simultaneous", 
                               verbose=FALSE, diagonal=TRUE, 
                               prior=list(pi=0.5), eps=1e-8){
  ##
  stop("linear not ready.")
  ## check prior
  if(!all(c("pi")%in%names(prior)))
    stop("Prior should be a list with elements (pi=numeric)")
  ##
  ## verbosity
  cat2 <- function(x)NULL
  if(verbose) cat2<-cat
  ##
  N <- length(y)
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
  }
  else{
    mb <- 0
    P <- 0
    Z <- Matrix(0)
  }
  ##
  ## Then the SS coef's:
  nongaussian <- !is.null(X)
  if(nongaussian){
    ## we use transposes
    X <- t(X)
    K <- nrow(X)
    D <- X%*%y   
    m0 <- rep(0, K)
    m <- rep(0, K)
    tau <- rep(1, K)
    pi <- rep(prior$pi, K)
    Ftilde <- rep(0, K)
  }
  else{
    K <- 0
    X <- Matrix(0)
    D <- Matrix(0)
    pi <- 0
  }
  cat2("Gaussian coefficients:", P, "\n")
  cat2("SS coefficients:", K, "\n")
  ##
  ## helper functions:
  lambda <- function(l)  (logit(l)-0.5)/(2*l)
  ## 
  ## initialize.
  Et  <- 1
  ##################################
  ## The optimization loop:
  loop <- TRUE
  while(loop){
    ## collect for convergence diagnostics:
    piold <- pi
    ##
    ## update SS parameters
    if(nongaussian){
      XtildeT <- la*t(X) 
      Btilde <- X %*% XtildeT
      ## Gaussian coefficient contribution to the mix
      if(gaussian){
        GtildeT <- Z %*% XtildeT
        Ftilde <- t(mb)%*%GtildeT
      }
      for(k in 1:K){
        ## slab
        tau[k] <- 2*Btilde[k,k] + Et
        Bk     <- t(pi[-k]*m[-k])%*%Btilde[-k,k]
        m[k]   <- (0.5*D[k] + Et*m0[k] - 2*Bk - 2*Ftilde[k] )/tau[k]
        ## spike
        q1 <- log(prior$pi) + 0.5*tau[k]*m[k]^2 - 0.5*Et*m0[k]^2
        q0 <- log(1-prior$pi) 
        u <- as.numeric(q1-q0)
        pi[k] <- logit(u)
      }
      ## expectation of coef and their squares. Note: diagonal cov assumption.
      w <- pi*m
      w2 <- pi*( 1/tau+m^2 )
      ## then we update the shared hyper-precision for the slabs:
      alfa <- 10 + K/2
      beta <- 10 + 0.5*sum(w2 - 2*w*m0 + m0^2)
      Et <- as.numeric(  alfa/beta  )
    }
    ##
    ## then update the gaussian coefficients, if they exist
    if(gaussian){
      ZtildeT <- la*t(Z) ## multiply each column i with la[i]
      Htilde <- Z%*% ZtildeT
      Sbi <- 2*Htilde + Sbi0
      ## this makes things a lot faster but is also biased
      if(diagonal) Sbi  <- Diagonal(x=diag(Sbi))
      Sb <- solve(Sbi)
      ## non gaussian contribution to the mix
      if(nongaussian) Etilde <- GtildeT%*%w
      else Etilde <- 0
      mb <- Sb%*% ( 0.5*C - 2*Etilde + Sbi0%*%mb0  )
      mb2 <- diag(Sb)+mb^2
    }
    ##
    ## then we update the logistic link variational parameters
    for(i in 1:N){
      xii <- 0
      if(nongaussian){
        Bi <- X[,i]%*%t(X[,i])
        xii <- xii + t(w)%*%Bi%*%(w) - t(w^2-w2)%*%diag(Bi)
      }
      if(gaussian){
        Oi <- Z[,i]%*%t(Z[,i])
        xii <- xii + t(mb)%*%Oi%*%(mb) - t(mb^2-mb2)%*%diag(Oi) 
        if(nongaussian){
          Yi <- X[,i]%*%t(Z[,i])
          xii <- xii + 2*t(w)%*%Yi%*%mb
        }
      }
      xi[i] <- xii[1]
    }
    xi<-sqrt(xi) 
    ##
    ## loop finishing...
    ####################################
    d1 <- abs(piold-pi)
    d2 <- abs(xiold-xi)
    loop <- (d<-max(d1, d2) )> eps
    cat2("d:", d, "           \r")
  }
  cat2("\n")
  ## done.
  w<-b<-NULL
  if(gaussian) b<-list(m=mb, S=Sb, Si=Sbi)
  if(nongaussian) w<-list(m=m, tau=tau, pi=pi, xi=xi, hyper=c(alpha=alfa, beta=beta))
  else b <- NULL
  list(w=w, b=b)
}
### End of VB linear regression SS function.



