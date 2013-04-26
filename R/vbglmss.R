#' Variational approximation of Spike and Slab GLM
#'
#' The function fits the GLM model in Bayesian framework with some covariate coefficients
#' having the spike and slab prior. 
#' @param formula an object of class \link{formula}. See Details.
#' @param family the likelihood family to be used in the model, see \link{family}.
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the variables in the model. 
#' If not found in data, the variables are taken from environment(formula), much
#' like in \link{glm}.
#' @param intercept should the intercept term be explicitely added.
#' @param y Alternative for formula&data: The response
#' @param X Alternative for formula&data: The SS covariates
#' @param Z Alternative for formula&data: The non-SS (regular Gaussian) covariates
#' @param method Use 'simultaneous' or 'sequential'. See Details.
#' @param verbose print additional runtime messages
#' @param eps the approximation loop convergence threshold
#' @details
#' The formula should use the special tags 'ss' for those covariates that are set 
#' the spike and slab prior. Without it, we assume regular Gaussian prior.
#' 
#' Using the y,X,Z variables is useful to avoid overhead in formula parsing when
#' data is large. This is then equal to the formula y ~ ss(X) + Z
#' @return 
#'  An object of class \link{vbglmss}.
#' @export

vbglmss <- function(formula, family=gaussian, data=NULL, intercept=FALSE,
                    y=NULL, X=NULL, Z=NULL, 
                    ...) {
  
  if(is.null(y) | (is.null(X) & is.null(Z)) )
    stop("formula not yet supported, use y, X parameters Z instead.")
  
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  
  fitter <- switch(family$family, 'binomial'=vbglmss.fit.logistic,
                                  'gaussian'=vbglmss.fit.linear)
  
  result <- fitter(y=y, X=X, Z=Z, ...)   
  
  result$family <- family
  result$call <- match.call()
  class(result) <- c("vbglmss")
  result
}

#' Print method for result object of \link{vbglmss}
#'
#' @param x the result object.
#' @S3method print vbglmss

print.vbglmss<-function(x, ...) {
  print(x)
}