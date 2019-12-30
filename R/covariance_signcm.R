#' signcm
#'
#' calculates the spatial sign covariance matrix (:= SCM)
#'
#' @param X : N x p matrix
#'
#' @param center : boolean. If True center data using spatial median
#'                Default = False
#'
#' @return C : Spatial covariance matrix
#' @return smed0 : spatial median. Computed only if center = True
#'
#' @examples
#'
#' x <- matrix(rnorm(9), 3, 3)
#'
#' signcm(x)
#' @export
signcm <- function(x, center = F){
  EPS <- 1e-6

  if(center){
    smed0 <- Rrobustsp::spatmed(x)
    x <- sweep(x, 1, smed0)
  }
  else {
    smed0 <- NULL
  }

  len <- sqrt(rowSums(x * Conj(x)))

  x[len != 0, ] <- x

  len <- len[len != 0]

  n <- nrow(x)

  len[len < EPS] <- EPS

  x <- sweep(x, 1, len, FUN = '/')

  C <- t(x) %*% Conj(x) / n

  return(list(C, smed0))
}
