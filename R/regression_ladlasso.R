#' ladlasso
#'
#' ladlasso computes the LAD-Lasso regression estimates for given complex-
#' or real-valued data.  If number of predictors, p, is larger than one,
#' then IRWLS algorithm is used, otherwise a weighted median algorithm
#' (N > 200) or elemental fits (N<200).
#'
#' @param y numeric response N x 1 vector (real/complex)
#' @param X  numeric feature  N x p matrix (real/complex)
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
#' @return b1 numeric, the regression coefficient vector of size N
#' @return integr, iter number of iterations
#' @export
#'
#' @examples
ladlasso <- function(y, X, lambda, intcpt = T, b0 = NULL, reltol = 1e-8, printitn = 0, iter_max = 2000){
  N <- nrow(X)
  p <- ncol(X)

  # make matrix sparse
  X <- Matrix(X)

  if(intcpt) X <- cbind(matrix(1, N, 1), X)

  if(is.null(b0)) b0 <- qr.solve(X, y) # ginv(X) %*% y

  iter <- NULL
  if(printitn > 0) sprintf('Computing the solution for lambda = %.3f\n',lambda)

  # The case of only one predictor
  if(p == 1){
    if(!intcpt) b1 <- wmed( rbind(y / X, 0), rbind(abs(X), lambda))
    if(!is.complex(y) & N < 200 & intcpt){
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


    for(iter in 1:iter_max){
      resid <- abs(y - X %*% b0)
      resid[resid < 1e-6] <- 1e-6

      Xstar <- sweep_sparse(X, 1, resid, fun = "/")


      b1 <- solve(Matrix::t(Xstar) %*% X, (Matrix::t(Xstar) %*% y)@x)

      crit <- norm(b1-b0, type = "2") / norm(b0, type = "2")

      if(printitn > 0 & iter %% printitn) sprintf('ladlasso: crit(%4d) = %.9f\n',iter,crit)
      if(crit < reltol && iter > 10) break
      b0 <- b1
    }}

  return( list(c(b1), iter))
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
