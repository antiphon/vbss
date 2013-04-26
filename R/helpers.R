#' helper functions
#' 
logit <- function(x)1/(1+exp(-x))

#' helper functions
#' 
lambda <- function(l)  (logit(l)-0.5)/(2*l)
