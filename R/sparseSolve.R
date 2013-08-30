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
  r<-try(v<-{
      L <- Cholesky(A, LDL=FALSE, super=!TRUE, perm=FALSE)
      solve(L, system="A")
      }, silent=TRUE)
  if("try-error" %in% is(r)) v <- solve(A)
  v
}