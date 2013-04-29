\name{vbglmss.mixedSS}
\alias{vbglmss.mixedSS}
\title{Variational approximation of Spike and Slab regression}
\usage{
  vbglmss.mixedSS(y, Xy, Zy, h, Xh, Zh,
    prior = list(theta = list(), beta = list(), binary = list(), gaussian = list()),
    verbose = FALSE, damp = 0, eps = 0.01, max.iter = 100)
}
\arguments{
  \item{y}{Binary observations}

  \item{Xy}{SS covariates for y}

  \item{Zy}{G covariates for y}

  \item{h}{Gaussian observations}

  \item{Xh}{SS covariates for h}

  \item{Zh}{G covariates for h}

  \item{prior}{A list of priors. See details.}

  \item{verbose}{print additional runtime messages}

  \item{damp}{Damping coefficient in SS mean update. damp *
  old + (1-damp)*new.}

  \item{eps}{the approximation loop convergence threshold}

  \item{max.iter}{--}
}
\value{
  A list with posteriors.
}
\description{
  The function fits a regression model where observations
  can be from logistic and/or gaussian likelihood. The
  coefficients are shared, and can have Spike and Slab (SS)
  and/or Gaussian (G) priors.
}
\details{
  TODO
}

