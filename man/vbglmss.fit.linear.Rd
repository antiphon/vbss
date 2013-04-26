\name{vbglmss.fit.linear}
\alias{vbglmss.fit.linear}
\title{Fit VB Spike and Slab linear regression}
\usage{
  vbglmss.fit.linear(y, X, Z, method = "simultaneous",
    verbose = FALSE, diagonal = TRUE,
    prior = list(pi = 0.5), eps = 1e-08)
}
\arguments{
  \item{y}{response.}

  \item{X}{the SS covariates}

  \item{Z}{the non-SS (regular Gaussian) covariates}

  \item{method}{Use 'simultaneous' or 'sequential'. See
  Details.}

  \item{verbose}{print additional runtime messages}

  \item{diagonal}{use diagonal covariance matrix for
  regular coefficients?}

  \item{eps}{the approximation loop convergence threshold}
}
\description{
  The workhorse for linear (Gaussian family)
  \link{vbglmss}.
}

