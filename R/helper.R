sweep_sparse <- function(x, margin, stats, fun = "-") {
  f <- match.fun(fun)
  if (margin == 1) {
    idx <- x@i + 1
  } else {
    idx <- x@j + 1
  }
  x@x <- f(x@x, stats[idx])
  return(x)
}

#' repmat
#'  matlab like repmat, vector v is replicate s times columnwise into a matrix
#'@export
repmat <- function(v, s) matrix(v, nrow=length(v), ncol=s)

#' mat_sign
#'
#' sign that also works for complex numbers
#' @export
mat_sign <- function(x){
  if(is.complex(x)) res <- x / abs(x) else res <- sign(x)
  return(res)
}

orth <- function(x){
  # rankMatrix is from Matrix library
  svd(x)$u[,1:rankMatrix(x)[[1]]]}

#' sqrtm
#'
#' computes square root of matrix via Jordan Normal form
#'
#' @export
sqrtm <- function(mat){
  X <- eigen(mat)
  T <- X$vectors
  vals <- sqrt(as.complex(X$values))
  J <- diag(vals)
  return(T %*% J %*% solve(T))
}

#' inf_norm
#'
#' norm(mat , type = 'I') discards
#' the imaginary part of complex matrices
#'
inf_norm <- function(mat){
  return(max(rowSums(abs(mat))))
}
