\name{vbglmss.fit.logistic}
\alias{vbglmss.fit.logistic}
\title{Fit VB Spike and Slab logistic regression}
\usage{
  vbglmss.fit.logistic(y, X = NULL, Z = NULL,
    prior = list(pi = 0.5, tau = list(a = 10, b = 10)),
    verbose = TRUE, ...)
}
\arguments{
  \item{y}{binary response.}

  \item{X}{the SS covariates.}

  \item{Z}{the non-SS (regular Gaussian) covariates.}

  \item{verbose}{print additional runtime messages}

  \item{eps}{the approximation loop convergence threshold}

  \item{...}{further parameters passed to sub-functions.
  This function makes the decision of which sub-functions
  need to be called based on the presense/absence of X and
  Z and parameter values.}
}
\description{
  The workhorse for logistic \link{vbglmss}.
}

