#' Variational approximation of Spike and Slab GLM.
#'
#' @param formula an object of class \link{formula}. See Details.
#' @param family the likelihood family to be used in the model, see \link{family}.
#' @param data an optional data frame, list or environment (or object coercible 
#' by as.data.frame to a data frame) containing the variables in the model. 
#' If not found in data, the variables are taken from environment(formula), much
#' like in \link{glm}.
#' @param intercept should the intercept term be explicitely included.
#' @param y Alternative for formula&data: The response
#' @param X Alternative for formula&data: The SS covariates
#' @param Z Alternative for formula&data: The non-SS (regular Gaussian) covariates
#' @param verbose print additional runtime messages
#' @param eps the approximation loop convergence threshold
#' @details
#' Jee joo
#' @return 
#'  An object of class \link{vbglmss}.
#' @export

vbglmss <- function(formula, family=gaussian, data=NULL, intercept=TRUE,
                    y=NULL, X=NULL, Z=NULL, 
                    verbose=FALSE, eps=1e-8
                    ) {
     
  
  
  
}

#' Print method for result object of \link{vbglmss}
#'
#' @param x the result object.
#' @S3method print vbglmss

print.vbglmss<-function(x, ...) {
  print(x)
}