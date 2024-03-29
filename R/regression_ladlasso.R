#' ladlasso
#'
#' ladlasso computes the LAD-Lasso regression estimates for given complex-
#' or real-valued data.  If number of predictors, p, is larger than one,
#' then IRWLS algorithm is used, otherwise a weighted median algorithm
#' (N > 200) or elemental fits (N<200).
#'
#' @param y numeric response N x 1 vector (real/complex)
#' @param X sparse matrix, numeric feature  N x p matrix (real/complex)
#' @param lambda numeric, non-negative penalty parameter
#' @param intcpt numeric optional initial start of the regression vector for
#'        IRWLS algorithm. If not given, we use LSE (when p>1).
#' @param b0 numeric optional initial start of the regression vector for
#'           IRWLS algorithm. If not given, we use LSE (when p>1).
#' @param reltol
#' Convergence threshold for IRWLS. Terminate when successive
#'           estimates differ in L2 norm by a rel. amount less than reltol.
#' @param printitn print iteration number (default = 0, no printing)
#' @param iter_max number of iterations \cr
#' default = 2000
#'
#' @return b1: numeric, the regression coefficient vector of size N
#' @return iter: integer, number of iterations
#'
#' @examples
#'
#' ladlasso(rnorm(8), matrix(rnorm(8*3)*5+1, 8, 3), 0.5)
#'
#' @note
#'
#' File location: regression_ladlasso.R
#' @export
ladlasso <- function(y, X, lambda, intcpt = T, b0 = NULL, reltol = 1e-8, printitn = 0, iter_max = 2000){
  N <- nrow(X)
  p <- if(is.null(ncol(X))) 1 else ncol(X)

  # make matrix sparse
  # X <- sparseMatrix(i = row(X)[row(X) != 0], j = col(X)[col(X) != 0], x=c(X))

  if(intcpt) X <- cbind(matrix(1, N, 1), X)

  if(is.null(b0)) b0 <- qr.solve(X, y) # ginv(X) %*% y

  iter <- NULL
  if(printitn > 0) sprintf('Computing the solution for lambda = %.3f\n',lambda)

  # The case of only one predictor
  if(p == 1){
    if(!intcpt){
      b1 <- wmed( rbind(y / X, 0), rbind(abs(X), lambda))
      return( list('b1' = c(b1), 'iter' = iter))}
    if(!is.complex(y) && N < 200 && intcpt){
      if(lambda == 0){
        b <- elemfits(X[,2], y) # b is a matrix
        b <- b[[1]]
      }else{
        b <- elemfits(c(X[,2], 0), c(y, lambda))
        b <- b[[1]]
      }}
    res <- colSums(abs(repmat(y, ncol(b)) - X %*% b))
    indx <- which.min(res)
    b1 <- b[,indx]
  }
  else {
    # use IRWLS always when p > 1
    if(printitn > 0) print('Starting the IRWLS algorithm..\n')
    if(lambda > 0){
      y <- c(y, rep(0, p))

      # slow
      if(intcpt) X <- rbind(X, cbind(rep(0, p), diag(lambda, p, p))) else X <- rbind(X, diag(lambda, p, p))
    }

    if(class(X)[1] == "matrix"){
      sweep2 <- function(x, y){sweep(x, MARGIN = 1, y, FUN = '/')}
      solve2 <- function(Xstar, X, y){solve(t(Xstar) %*% X, (t(Xstar) %*% y))}
      }
    else{
      sweep2 <- function(x, y){sweep_sparse(x, margin = 1, y, fun = '/')}
      solve2 <- function(Xstar, X, y){solve(Matrix::t(Xstar) %*% X, (Matrix::t(Xstar) %*% y)@x)}
    }

    for(iter in 1:iter_max){
      resid <- abs(y - X %*% b0)
      resid[resid < 1e-6] <- 1e-6

      # make if else for dense matrices
      Xstar <- sweep2(X, resid)

      b1 <- solve2(Xstar, X, y)

      crit <- norm(b1-b0, type = "2") / norm(b0, type = "2")

      if(printitn > 0 & iter %% printitn) sprintf('ladlasso: crit(%4d) = %.9f\n',iter,crit)
      if(crit < reltol && iter > 10) break
      b0 <- b1
    }}

  return( list('b1' = c(b1), 'iter' = iter))
}

# fRcpp ----

# library(Rcpp)
#
# cppFunction(depends='RcppArmadillo', code='
#             arma::mat fRcpp (arma::mat Xstar, arma::mat X, arma::mat y) {
#             arma::mat betahat ;
#             betahat = (Xstar.t() * X ).i() * (Xstar.t() * y) ;
#             return(betahat) ;
#             }
#             ')
