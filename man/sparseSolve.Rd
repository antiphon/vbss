\name{sparseSolve}
\alias{sparseSolve}
\title{Invert matrix}
\usage{
  sparseSolve(A)
}
\arguments{
  \item{A}{square matrix to invert.}
}
\value{
  Inverse of A.
}
\description{
  Invert matrix
}
\details{
  We use Cholesky decomposition for inversion, so A should
  be symmetric.
}

