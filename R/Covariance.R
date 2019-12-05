#' spatmed
#'
#' Computes the spatial median based on (real or complex) data matrix X.
#'
#' @param X : numeric, matrix of size N x p
#' @param printitn : 0 or 1. print iteration number. Default = 0.
#'
#' @return smed : spatial median estimate
#' @export
spatmed <- function(X, printitn = 0){
  len <- rowSums(X * Conj(X))
  X <- X[len != 0,]
  n <- nrow(X)

  if(is.complex(X)) smed0 <- colMeans(X) else smed0 <- median(X)

  norm0 <- norm(smed0, type = "2")
  iter_max <- 500
  eps <- 1e-6
  tol <- 1e-5

  for(iter in 1:iter_max){
    Xc <- X - smed0

    len <- sqrt(rowSums(Xc * Conj(Xc)))
    len[Re(len) < eps] <- eps

    Xpsi <- Xc / len

    update <- colSums(Xpsi) / sum(1 / len)

    smed <- smed0 + update

    dis <- norm(update, type = "2") / norm0

    if(printitn >0 & iter %% printitn) sprintf('At iter = %3d, dis=%.7f\n',iter,dis)

    if(dis <= tol) break

    smed0 <- smed

    norm0 <- norm(smed, type = "2")
  }

  return(smed)
}
