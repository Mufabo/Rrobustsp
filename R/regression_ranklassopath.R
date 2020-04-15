#' ranklassopath
#'
#' ranklassopath computes the rank LAD-Lasso regularization path (over grid
#'                                                                  of penalty parameter values). Uses IRWLS algorithm.
#'
#' @param         y: Numeric data vector of size N (output, respones)
#' @param       X: Numeric data matrix of size N x p. Each row represents one
#'           observation, and each column represents one predictor (feature).
#' @param       L: Positive integer, the number of lambda values on the grid to be
#'           used. The default is L=120.
#' @param     eps: Positive scalar, the ratio of the smallest to the
#'           largest Lambda value in the grid. Default is eps = 10^-3.
#' @param  reltol: Convergence threshold for IRWLS. Terminate when successive
#'          estimates differ in L2 norm by a rel. amount less than reltol.
#' @param printitn: print iteration number (default = F, no printing)
#'
#' @return          B: Fitted RLAD-Lasso regression coefficients, a p-by-(L+1) matrix,
#'           where p is the number of predictors (columns) in X, and L is
#'           the  number of Lambda values.
#' @return      B0: estimates values of intercepts
#' @return   stats: structure with following fields:
#'             Lambda = lambda parameters in ascending order
#'             GMeAD = Mean Absolute Deviation (MeAD) of the residuals
#'             gBIC = generalized Bayesian information criterion (gBIC) value
#'                  for each lambda parameter on the grid.
#'
#'@examples
#' data('prostate')
#'
#' X <- prostate$X
#' y <- c(prostate$y)
#'
#' namess <- unlist(prostate$names)
#'
#' n <- nrow(X)
#' p <- ncol(X)
#'
#' Xone <- cbind(rep(1,n), X)
#' LSE <- qr.solve(Xone, y) # Least squares estimate
#'
#' GRlen <- 120
#'
#' yout <- y
#' yout[1] <- yout[1] + 55
#' ranklassopath(yout, X)
#' @note
#' File location : regression_ranklassopath.R
#' @export
ranklassopath <- function(y, X, L = 120, eps = 1e-3, reltol = 1e-7, printitn = F){
  n <- nrow(X)
  p <- ncol(X)
  intcpt <- F

  B <- repmat(1:n, n)

  A <- t(B)
  a <- A[A < B]
  b <- B[A < B]

  Xtilde <- X[a,] - X[b,]
  ytilde <- y[a] - y[b]

  lam0 <- norm(t(Xtilde) %*% mat_sign(ytilde), type = "I")

  lamgrid <- eps^((0:L)/L) * lam0

  B <- diag(0, p, L + 1)
  B0 <- numeric(L + 1)
  b_init <- numeric(p)

  for(jj in 1:(L + 1)){
    B[, jj] <- ladlasso(ytilde, Xtilde, lamgrid[jj], intcpt, b_init, reltol, printitn)[[1]]
    b_init <- B[, jj]
    r <- y - X %*% b_init
    if(is.complex(X)) B0[jj] <- spatmed((r[a] + r[b]) / 2) else B0[jj] <- median((r[a] + r[b]) / 2)

  }

  B[abs(B) < 1e-7] <- 0

  DF <- sum(abs(B) != 0)

  # Compute the generalized BIC (gBIC) values
  Rmat <- repmat(ytilde, L + 1) - Xtilde %*% B # matrix of residuals
  N <- (n - 1) / 2
  GMeAD <- (sqrt(pi) / 2) * colMeans(abs(Rmat)) # Gini's dispersion
  GmeAD <- GMeAD * sqrt(n / (n - DF -1))

  gBIC <- 2 * n * log(GMeAD) + DF * log(n)

  Lambda <- lamgrid

  stats <- list('DF' = DF, 'GMeAD' = GMeAD, 'gBIC' = gBIC, 'Lambda' = Lambda)

  return(list(B, B0, stats))
}
