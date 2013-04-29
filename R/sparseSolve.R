#' Invert matrix
#'
#' @param A square matrix to invert.
#' @details
#' We use Cholesky decomposition for inversion, so A should be symmetric.
#' @return 
#'  Inverse of A.
#' @import Matrix
#' @export

sparseSolve <- function(A) {
  #L <- Cholesky(A, LDL=FALSE, super=TRUE, perm=FALSE)
  #solve(L, system="A")
  solve(A)
}