\name{vbglmss}
\alias{vbglmss}
\title{Variational approximation of Spike and Slab GLM.}
\usage{
  vbglmss(formula, family = gaussian, data = NULL,
    intercept = TRUE, y = NULL, X = NULL, Z = NULL,
    verbose = FALSE, eps = 1e-08)
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
  included.}

  \item{y}{Alternative for formula&data: The response}

  \item{X}{Alternative for formula&data: The SS covariates}

  \item{Z}{Alternative for formula&data: The non-SS
  (regular Gaussian) covariates}

  \item{verbose}{print additional runtime messages}

  \item{eps}{the approximation loop convergence threshold}
}
\value{
  An object of class \link{vbglmss}.
}
\description{
  Variational approximation of Spike and Slab GLM.
}
\details{
  Jee joo
}

