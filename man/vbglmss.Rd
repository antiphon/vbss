\name{vbglmss}
\alias{vbglmss}
\title{Variational approximation of Spike and Slab GLM}
\usage{
  vbglmss(formula, family = gaussian, data = NULL,
    intercept = FALSE, y = NULL, X = NULL, Z = NULL, ...)
}
\arguments{
  \item{formula}{an object of class \link{formula}. See
  Details.}

  \item{family}{the likelihood family to be used in the
  model, see \link{family}.}

  \item{data}{an optional data frame, list or environment
  (or object coercible by as.data.frame to a data frame)
  containing the variables in the model. If not found in
  data, the variables are taken from environment(formula),
  much like in \link{glm}.}

  \item{intercept}{should the intercept term be explicitely
  added.}

  \item{y}{Alternative for formula&data: The response}

  \item{X}{Alternative for formula&data: The SS covariates}

  \item{Z}{Alternative for formula&data: The non-SS
  (regular Gaussian) covariates}

  \item{method}{Use 'simultaneous' or 'sequential'. See
  Details.}

  \item{verbose}{print additional runtime messages}

  \item{eps}{the approximation loop convergence threshold}
}
\value{
  An object of class \link{vbglmss}.
}
\description{
  The function fits the GLM model in Bayesian framework
  with some covariate coefficients having the spike and
  slab prior.
}
\details{
  The formula should use the special tags 'ss' for those
  covariates that are set the spike and slab prior. Without
  it, we assume regular Gaussian prior.

  Using the y,X,Z variables is useful to avoid overhead in
  formula parsing when data is large. This is then equal to
  the formula y ~ ss(X) + Z
}

