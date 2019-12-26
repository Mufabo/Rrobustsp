# hublassopath ----

#'hublassopath
#'
#'hublassopath computes the M-Lasso regularization path (over grid
#'of penalty parameter values) using Huber's loss function.
#'
#'@param          y: Numeric data vector of size N x 1 (output, respones)
#'@param          X: Numeric data matrix of size N x p. Each row represents one
#'             observation, and each column represents one predictor (feature)
#'             columns are standardized to unit length.
#'@param          c: Threshold constant of Huber's loss function (optional;
#'                                                                       otherwise use default value)
#'@param       intcpt: Logical (true/false) flag to indicate if intercept is in the
#'             regression mode. Default is true.
#'@param        eps: Positive scalar, the ratio of the smallest to the
#'             largest Lambda value in the grid. Default is eps = 10^-3.
#'@param         L : Positive integer, the number of lambda values EN/Lasso uses.
#'             Default is L=120.
#'@param     reltol : Convergence threshold for IRWLS. Terminate when successive
#'             estimates differ in L2 norm by a rel. amount less than reltol. default: 1e-05
#'@param   printitn: print iteration number (default = 0, no printing)
#'
#'@return       B    : Fitted M-Lasso regression coefficients, a p-by-(L+1) matrix,
#'            where p is the number of predictors (columns) in X, and L is
#'            the  number of Lambda values.
#'@return       B0 : estimates values of intercepts
#'@return   stats  : structure with following fields:
#'\itemize{
#'               \item Lambda = lambda parameters in ascending order
#'             \item sigma = estimates of the scale (a (L+1) x 1 vector)
#'             \item gBIC = generalized Bayesian information criterion (gBIC) value
#'                   for each lambda parameter on the grid.}
#'@export
hublassopath <- function(y, X, c = NULL, L = 120, eps =10^-3, intcpt = T, reltol = 1e-5, printitn = 0){
  n <- nrow(X)
  p <- ncol(X)

  if(is.complex(y)) real_data <- F else real_data <- T

  # Default: approx 95 efficiency for Gaussian errors
  if(is.null(c)){if(real_data) c <- 1.3414 else c <- 1.215}

  tmp <- hubreg(y, matrix(1, n, 1), c)
  locy <- tmp[[1]]
  sig0 <- tmp[[2]]

  if(intcpt){
    # center data
    ync <- y
    Xnc <- X

    meanX <- colMeans(X)
    X <- sweep(X, 2, meanX) #X - meanX
    y <- y - locy
  }

  # standardize the predictor to unit norm columns
  sdX <- sqrt(colSums(X * Conj(X)))
  X <- sweep(X, 2, sdX, FUN = '/')

  # compute the smallest penalty value yielding a zero solution
  yc <- psihub(y / sig0, c) * sig0
  lam0 <- norm(t(X) %*% yc, type = "i")

  lamgrid <- eps^((0:L) / L) * lam0
  B <- matrix(0, p, L+1)
  sig <- matrix(0, 1, L+1)
  sig[1] <- sig0

  for(jj in 1:L){
    tmp <- hublasso(y, X, c, lamgrid[jj+1], B[,jj], sig[jj], reltol, printitn)
    B[, jj+1] <- tmp[[1]]
    sig[jj+1] <- tmp[[2]]
  }

  B[abs(B) < 5e-8] <- 0
  DF <- colSums(abs(B) != 0)
  con <- sqrt((n/(n-DF-1)))

  if(n > p) gBIC <- 2 * n * log(sig * con) + DF * log(n) else gBIC <- NULL

  B <- sweep(B, 1, sdX, FUN = '/')

  # compute the intercept if in the model
  if(intcpt) B0 <- locy - c(meanX %*% B) else B <- NULL

  stats <- list(gBIC, sig, lamgrid)
  names(stats) <- c('gBIC', 'sigma', 'Lambda')

  return( list('B' = B, 'B0' = B0, 'stats' = stats))
}
