#'   [beta, w] = elemfits(x,y)
#'   elemfits compute the Nx(N-1)/2 elemental fits, i.e., intercepts b_{0,ij}
#'   and slopes b_{1,ij}, that define a line y = b_0+b_1 x that passes through
#'   the data points (x_i,y_i) and (x_j,y_j), i<j, where i, j in {1, ..., N}.
#'   and the respective weights | x_i - x_j |
#'
#'    @param   y : (numeric) N x 1 vector of real-valued outputs (response vector)
#'    @param   x : (numeric) N x 1 vector of inputs (feature vector)
#'
#'    @return    beta: (numeric) N*(N-1)/2 matrix of elemental fits
#'    @return   w: (numeric) N*(N-1)/2 matrix of weights
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




#'   enet computes the elastic net estimator using the cyclic co-ordinate
#'   descent (CCD) algorithm.
#'
#'     @param       y : (numeric) data vector of size N x 1 (output, respones)
#'           if the intercept is in the model, then y needs to be centered.
#'     @param    X : (numeric) data matrix of size N x p (input, features)
#'             Columns are assumed to be standardized (i.e., norm(X(:,j))=1)
#'             as well as centered (if intercept is in the model).
#'   @param   beta : (numeric) regression vector for initial start for CCD algorithm
#'   @param lambda : (numeric) a postive penalty parameter value
#'   @param alpha  : (numeric) elastic net tuning parameter in the range [0,1]. If
#'             not given then use alpha = 1 (Lasso)
#'  @param printitn: print iteration number (default = 0, no printing)
#'
#'   @return     beta    : (numberic) the regression coefficient vector
#'  @return iter  : (numeric) # of iterations
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

    if(is.nan(crit) | crit < 1e-4) {      break }

    betaold <- beta
    normb0 <- normb
  }

  return (list(beta, iter))
}



#'  [B, lamgrid, BIC, MSE] = enetpath(y, X,...)
#'  enethpath computes the elastic net (EN) regularization path (over grid
#' of penalty parameter values). Uses pathwise CCD algorithm.
#'
#'  @param       y : Numeric data vector of size N x 1 (output, respones)
#'  @param     X : Nnumeric data matrix of size N x p. Each row represents one
#'          observation, and each column represents one predictor (feature).
#'  @param intcpt: Logical flag to indicate if intercept is in the model
#'  @param alpha  : Numeric scalar, elastic net tuning parameter in the range [0,1].
#'            If not given then use alpha = 1 (Lasso)
#'   @param    eps: Positive scalar, the ratio of the smallest to the
#'           largest Lambda value in the grid. Default is eps = 10^-4.
#'   @param     L : Positive integer, the number of lambda values EN/Lasso uses.
#'            Default is L=100.
#'  @param printitn: print iteration number (default = 0, no printing)
#'
#'  @return    B    : Fitted EN/Lasso regression coefficients, a p-by-(L+1) matrix,
#'           where p is the number of predictors (columns) in X, and L is
#'           the  number of Lambda values. If intercept is in the model, then
#'           B is (p+1)-by-(L+1) matrix, with first element the intercept.
#'  @return stats  : structure with following fields:
#'              Lambda = lambda parameters in ascending order
#'            MSE = Mean squared error (MSE)
#'            BIC = Bayesian information criterion values for each lambda
#' @export
enetpath <- function(y, X, alpha=1,  L=100, eps=1e-4, intcpt=TRUE, printitn=0){
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
    B[,jj+1] <- enet(y, X, B[,jj], lamgrid[jj+1], alpha, printitn)[[1]]
  }

  B[abs(B)<5e-8] <- 0

  DF = colSums(abs(B) != 0)

  if(n>p){
      MSE <- colSums(abs(y - X %*% B)^2) * (1 / (n - DF - 1))
      BIC <- n * log(MSE) + DF * log(n)
  } else
  {
    MSE <- NULL
    BIC <- NULL
  }


  B <- B / sdX

  if(intcpt) B <- rbind(meany - meanX %*% B, B)

  stats <- matrix(c(MSE, BIC, lamgrid), nrow = L+1, ncol = 3, dimnames = list(NULL ,c('MSE', 'BIC', 'Lambda')))

  return (list(B, stats))
}



#'   hublasso computes the M-Lasso estimate for a given penalty parameter
#'   using Huber's loss function
#'
#'   @param       y: Numeric data vector of size N x 1 (output,respones)
#'   @param       X: Numeric data matrix of size N x p (inputs,predictors,features).
#'             Each row represents one observation, and each column represents
#'             one predictor
#'
#'   @param  lambda: positive penalty parameter value
#'   @param      b0: numeric initial start of the regression vector
#'   @param    sig0: numeric positive scalar, initial scale estimate.
#'   @param  c: Threshold constant of Huber's loss function
#'   @param  reltol: Convergence threshold. Terminate when successive
#'             estimates differ in L2 norm by a rel. amount less than reltol.
#'             Default is 1.0e-5
#'   @param  iter_max: int, default = 500. maximum number of iterations
#'   @param  printitn: print iteration number (default = 0, no printing)
#'
#'   @return        b0: regression coefficient vector estimate
#'   @return     sig0: estimate of the scale
#'   @return   psires: pseudoresiduals
#' @export
hublasso <- function(y, X, c = NULL, lambda, b0,  sig0, reltol = 1e-5, printitn = 0, iter_max = 500){
  n <- nrow(X)
  p <- ncol(X)

  complex_data <- is.complex(X)

  if(is.null(c)){if(complex_data) c <- 1.1214 else c <- 1.345}

  csq <- c^2

  # consistency factor for scale
  if(complex_data){
    qn <- pchisq(2 * csq, 2)
    alpha <- pchisq(2 * csq, 4) + csq * (1 - qn)
  } else{
    qn <- pchisq(csq, 1)
    alpha <- pchisq(csq, 3) + csq * (1 - qn)
  }

  con <- sqrt(n * alpha)
  betaold <- b0
  normb0 <- norm(matrix(b0), type = "2")

  for(iter in 1:iter_max){
    r <- y - X %*% b0
    psires <- psihub(r/sig0, c) * sig0
    # scale update
    sig1 <- norm(matrix(psires), type = "2") / con

    crit2 <- abs(sig0 - sig1)
    for(jj in 1:p){
      # update pseudoresidual
      psires <- psihub(r/sig1, c) * sig1
      b0[jj] <- soft_thresh(b0[jj] + X[,jj] %*% psires, lambda)
      r <- r + X[,jj] * (betaold[jj] - b0[jj])
    }

    normb <- norm(matrix(b0), type = "2")

    crit <- sqrt(normb0^2 + normb^2 - 2*Re(betaold %*% b0)) / normb0

    if(printitn > 0){
      r <- (y - X %*% b0) / sig1
      objnew <- sig1 * colSums(rhofun(r, c)) + (beta / 2) * n * sig1 + lambda * sum(abs(b0))
      sprinf('Iter %3d obj =  %10.5f crit = %10.5f crit2 = %10.5f\n',
             iter,objnew,crit,crit2)
    }

    if(is.na(crit) | crit < reltol){ break }

    sig0 <- sig1
    betaold <- b0
    normb0 <- normb
  }

  if(printitn > 0){
    # Check the M-Lasso estimating equations
    b0[abs(b0) < 5e-9] <- 0
    r <- y - X * b0
    s <- sign(b0)
    ind <- 1:p
    ind2 <- ind[s==0]
    psires <- psihub(r / sig1, c) * sig1
    s[ind2] <- X[, ind2] * psires / lambda
    # FP equation equal to zero
    fpeq <- -X * psires + lambda * s
    sprintf('lam = %.4f it = %d norm(FPeq1)= %.12f abs(FPeq2)=%.12f\n',
            lambda,iter, norm(fpeq),sig1-norm(psires)/con)
  }

  return (list(b0, sig0, psires))
}



#'   hublassopath computes the M-Lasso regularization path (over grid
#'                                                            of penalty parameter values) using Huber's loss function
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

  # tmp <- hubreg(y, rep(1, n), c)
  # TODO continue when done with hubreg
}



#'   hubreg computes the joint M-estimates of regression and scale using
#'   Huber's criterion. Function works for both real- and complex-valued data.
#'
#'
#'@param   y: Numeric data vector of size N x 1 (output, respones)
#'@param   X: Numeric data matrix of size N x p. Each row represents one
#'             observation, and each column represents one predictor (feature).
#'             If the model has an intercept, then first column needs to be a
#'             vector of ones.
#'@param    c: numeric threshold constant of Huber's function
#'@param    sig0: (numeric) initial estimator of scale [default: SQRT(1/(n-p)*RSS)]
#'@param    b0: initial estimator of regression (default: LSE)
#'@param    printitn: print iteration number (default = 0, no printing)
#'@param    iter_max: maximum number of iterations. \cr default = 2000
#'@param    errortol: ERROR TOLERANCE FOR HALTING CRITERION. \cr default = 1e-5
#'
#'@return    b1: the regression coefficient vector estimate
#'@return    sig1: the estimate of scale
#'@return    iter: the # of iterations
#'@export
hubreg <- function(y, X, c = NULL, sig0 = NULL, b0 = NULL, printitn = 0, iter_max = 2000, errortol = 1e-5){
  n <- nrow(X)
  p <- ncol(X)

  if(is.complex(y)) real_data <- F else real_data <- T
  # Default: approx 95 efficiency for Gaussian errors
  if(is.null(c)){if(real_data) c <- 1.3415 else c <- 1.215}

  if(is.null(b0)) b0 <- solve(X, y)

  # use unbiased residual sum of squares as initial estimate of scale
  if(is.null(sig0)){sig0 <- norm(matrix(y - X %*% b0), type = "2") /
    sqrt(n - p)}

  csq <- c^2

  # consistency factor for scale
  if(complex_data){
    qn <- pchisq(2 * csq, 2)
    alpha <- pchisq(2 * csq, 4) + csq * (1 - qn)
  } else{
    qn <- pchisq(csq, 1)
    alpha <- pchisq(csq, 3) + csq * (1 - qn)
  }

  Z <- pinv(X)

  con <- sqrt((n - p) * alpha)

  for(iter in 1:iter_max){
    # Step 1: Update residual
    r <- y - X %*% b0
    psires <- psihub(r / sig0, c) * sig0

    # Step 2: Update the scale
    sig1 <- norm(matrix(psires), type ="2") / con

    # Step 3: Update the pseudo-residual
    psires <- psihub(r / sig1, c) * sig1

    # Step 4: regress X on pseudo-residual
    update <- Z * psires

    # Step 5: update the beta
    b1 <- b0 + update
  }
}
