# % [beta, w] = elemfits(x,y)
# % elemfits compute the Nx(N-1)/2 elemental fits, i.e., intercepts b_{0,ij}
# % and slopes b_{1,ij}, that define a line y = b_0+b_1 x that passes through
# % the data points (x_i,y_i) and (x_j,y_j), i<j, where i, j in {1, ..., N}.
# % and the respective weights | x_i - x_j |
#   % INPUT:
#   %    y : (numeric) N x 1 vector of real-valued outputs (response vector)
# %    x : (numeric) N x 1 vector of inputs (feature vector)
# % OUTPUT:
#   %  beta: (numeric) N*(N-1)/2 matrix of elemental fits
# %     w: (numeric) N*(N-1)/2 matrix of weights
#' @export
elemfits <- function(y, x){
  N <- length(x)
  B <- matrix(1:N, N, N)
  a <- A[A<B]
  b <- B[A<B]
  w <- x[a] - x[b]

  beta <- matrix(0, 2, length(a))

  beta[2,] <- (y[a] - y[b])/w
  beta[1,] <- (x[a] * y[b] - x[b]*y[a])/w
  return (list(beta, abs(w)))
}

#' % enet computes the elastic net estimator using the cyclic co-ordinate
#' % descent (CCD) algorithm.
#' %
#' % INPUT:
#'   %       y : (numeric) data vector of size N x 1 (output, respones)
#' %         if the intercept is in the model, then y needs to be centered.
#' %       X : (numeric) data matrix of size N x p (input, features)
#' %           Columns are assumed to be standardized (i.e., norm(X(:,j))=1)
#' %           as well as centered (if intercept is in the model).
#' %    beta : (numeric) regression vector for initial start for CCD algorithm
#' %  lambda : (numeric) a postive penalty parameter value
#' %  alpha  : (numeric) elastic net tuning parameter in the range [0,1]. If
#' %           not given then use alpha = 1 (Lasso)
#' % printitn: print iteration number (default = 0, no printing)
#' % OUTPUT:
#'   %   beta    : (numberic) the regression coefficient vector
#' %   iter  : (numeric) # of iterations
#' @export
enet <- function(y, X, beta, lambda, alpha=1, printitn=0, itermax = 1000){
  # check for valid arguments

  p <- ncol(X)

  betaold <- beta
  normb0 <- norm(matrix(beta), type = "2")
  r <- y - X %*% beta

  lam1 <- alpha * lambda
  lam2 <- lambda * (1-alpha)

  const <- 1/(1+lam2)

  if(printitn > 0){
    cat('enet : using penalty lambda = ', round(lambda, digits=5))
  }

  for(iter in 1:itermax){
    for(jj in 1:p){
      beta[jj] <- const*soft_thresh(beta[jj]+ t(X[,jj]) %*% r, lam1)
      r <- r + X[,jj]*(betaold[jj]-beta[jj])
    }

    normb <- norm(matrix(beta), type = "2")

    crit <- sqrt(normb0^2 + normb^2 - 2 * Re(t(betaold) %*% beta))/normb

    if(!is.na(iter%%printitn) && (iter %% printitn) == 0){
      sprintf('enet: %4d  crit = %.8f\n',iter,crit)
    }

    if(crit < 1e-4) {      break }

    betaold <- beta
    normb0 <- normb
  }

  return (list(beta, iter))
}

#'  [B, lamgrid, BIC, MSE] = enetpath(y, X,...)
#'  enethpath computes the elastic net (EN) regularization path (over grid
#' of penalty parameter values). Uses pathwise CCD algorithm.
#'  INPUT:
#'        y : Numeric data vector of size N x 1 (output, respones)
#'       X : Nnumeric data matrix of size N x p. Each row represents one
#'          observation, and each column represents one predictor (feature).
#'   intcpt: Logical flag to indicate if intercept is in the model
#'   alpha  : Numeric scalar, elastic net tuning parameter in the range [0,1].
#'            If not given then use alpha = 1 (Lasso)
#'       eps: Positive scalar, the ratio of the smallest to the
#'           largest Lambda value in the grid. Default is eps = 10^-4.
#'        L : Positive integer, the number of lambda values EN/Lasso uses.
#'            Default is L=100.
#'  printitn: print iteration number (default = 0, no printing)
#'  OUTPUT:
#'      B    : Fitted EN/Lasso regression coefficients, a p-by-(L+1) matrix,
#'           where p is the number of predictors (columns) in X, and L is
#'           the  number of Lambda values. If intercept is in the model, then
#'           B is (p+1)-by-(L+1) matrix, with first element the intercept.
#'   stats  : structure with following fields:
#'              Lambda = lambda parameters in ascending order
#'            MSE = Mean squared error (MSE)
#'            BIC = Bayesian information criterion values for each lambda
#' @export
enetpath <- function(y, X, intcpt=TRUE, alpha=1, eps=1e-4, L=100, printitn=0){
  # TODO check for valid arguments

  n <- nrow(X)
  p <- ncol(X)

  # center the data if intcpt is True
  if(intcpt){
    meanX <- colMeans(X)
    meany <- mean(y)
    X <- X - meanX
    y <- y - meany
  }

  if (printitn>0) sprintf('enetpath: using alpha = %.1f \n', alpha)

  sdX = sqrt(colSums(X * Conj(X)))
  X <- X / sdX

  # smallest penalty value giving zero solution
  lam0 <- norm(t(X) %*% y , type="I")/alpha

  # grid of penalty values
  lamgrid <- eps^((0:L)/L) * lam0

  B <- matrix(0,p,L+1)

  for(jj in 1:L){
    B[,jj+1] <- enet(y, X, B[,jj], lamgrid[jj+1], alpha, printitn)
  }

  B[abs(B)<5e-8] <- 0

  DF = colSums(abs(B) != 0)

  if(n>p){
      MSE <- colSums(y - abs(X %*% B)^2) * (1 / (n - DF - 1))
      BIC <- n * log(MSE) + DF * log(n)
  } else
  {
    MSE <- NULL
    BIC <- NULL
  }


  B <- B / sdX
  if(intcpt) B <- cbind(meany - meanX %*% B, B)

  stats <- c(MSE, BIC, lamgrid)
  names(stats) <- c('MSE', 'BIC', 'Lambda')

  return (list(B, stats))
}
