# elemfits ----

#'elemental fits
#'
#'elemfits compute the Nx(N-1)/2 elemental fits, i.e., intercepts b_{0,ij}
#'and slopes b_{1,ij}, that define a line y = b_0+b_1 x that passes through
#'the data points (x_i,y_i) and (x_j,y_j), i<j, where i, j in {1, ..., N}.
#'and the respective weights | x_i - x_j |
#'
#'@param   y : (numeric) N  vector of real-valued outputs (response vector)
#'@param   x : (numeric) N  vector of inputs (feature vector)
#'
#'@return    beta: (numeric) N*(N-1)/2 matrix of elemental fits
#'@return   w: (numeric) N*(N-1)/2 matrix of weights
#'
#'@examples
#'
#'y <- c(0.5377   , 1.8339   ,-2.2588 ,   0.8622,    0.3188)
#'x <- c(-1.3077  , -0.4336,    0.3426   , 3.5784,    2.7694)
#'
#'elemfits(y, x)
#'
#'@note
#'File located in Regression.R
#'
#'@export
elemfits <- function(y, x){
  N <- length(x)
  B <- matrix(1:N, N, N)
  A <- t(B)
  a <- A[A<B]
  b <- B[A<B]
  w <- x[a] - x[b]

  beta <- matrix(0, 2, length(a))

  beta[2,] <- (y[a] - y[b])/w
  beta[1,] <- (x[a] * y[b] - x[b]*y[a])/w
  return (list(beta, abs(w)))
}


# enetpath ----

#' enetpath
#'
#' enethpath computes the elastic net (EN) regularization path (over grid
#' of penalty parameter values). Uses pathwise CCD algorithm.
#'
#'@param       y : Numeric data vector of size N x 1 (output, respones)
#'@param     X : Nnumeric data matrix of size N x p. Each row represents one
#'          observation, and each column represents one predictor (feature).
#'@param intcpt: Logical flag to indicate if intercept is in the model.
#'               Default = True
#'@param alpha  : Numeric scalar, elastic net tuning parameter in the range [0,1].
#'            If not given then use alpha = 1 (Lasso)
#'@param    eps: Positive scalar, the ratio of the smallest to the
#'           largest Lambda value in the grid. Default is eps = 10^-4.
#'@param     L : Positive integer, the number of lambda values EN/Lasso uses.
#'            Default is L=100.
#'@param printitn: print iteration number (default = 0, no printing)
#'
#'@return    B    : Fitted EN/Lasso regression coefficients, a p-by-(L+1) matrix,
#'           where p is the number of predictors (columns) in X, and L is
#'           the  number of Lambda values. If intercept is in the model, then
#'           B is (p+1)-by-(L+1) matrix, with first element the intercept.
#'@return stats  : list with following fields:
#'              \item{Lambda}{lambda parameters in ascending order}
#'            \item{MSE}{Mean squared error (MSE)}
#'            \item{BIC}{Bayesian information criterion values for each lambda}
#'@examples
#'y <- c(0.5377   , 1.8339   ,-2.2588 ,   0.8622,    0.3188)
#'x <- c(-1.3077  , -0.4336,    0.3426   , 3.5784,    2.7694)
#'
#'enetpath(y, matrix(c(x, x), nrow = 5, ncol = 2))
#'
#'@note
#'file in Regression.R
#'@export
enetpath <- function(y, X, alpha=1,  L=100, eps=1e-4, intcpt=TRUE, printitn=0){
  # TODO check for valid arguments

  n <- nrow(X)
  p <- ncol(X)

  # center the data if intcpt is True
  if(intcpt){
    meanX <- colMeans(X)
    meany <- mean(y)
    X <- sweep(X, 2, meanX) #X - meanX
    y <- y - meany
  }

  if (printitn>0) sprintf('enetpath: using alpha = %.1f \n', alpha)

  sdX = sqrt(colSums(X * Conj(X)))
  X <- sweep(X, 2, sdX, FUN = '/')

  # smallest penalty value giving zero solution
  lam0 <- norm(t(X) %*% y , type="I")/alpha

  # grid of penalty values
  lamgrid <- eps^((0:L)/L) * lam0

  B <- matrix(0, p, L+1)

  for(jj in 1:L){
    B[,jj+1] <- enet(y, X, B[,jj], lamgrid[jj+1], alpha, printitn)[[1]]
  }

  B[abs(B)<5e-8] <- 0

  DF = colSums(abs(B) != 0)

  if(n > p){
      MSE <- colSums(abs(y - X %*% B)^2) * (1 / (n - DF - 1))
      BIC <- n * log(MSE) + DF * log(n)
  } else
  {
    MSE <- NULL
    BIC <- NULL
  }


  B <- B / sdX

  if(intcpt) B <- rbind(meany - meanX %*% B, B)

  stats <- list(MSE, BIC, lamgrid)
  names(stats) <- c('MSE', 'BIC', 'Lambda')


  return (list(B, stats))
}


# hublasso ----

#'Lasso with Huber's loss function
#'
#'hublasso computes the M-Lasso estimate for a given penalty parameter
#'using Huber's loss function
#'
#'@param       y: Numeric data vector of size N x 1 (output,respones)
#'@param       X: Numeric data matrix of size N x p (inputs,predictors,features).
#'             Each row represents one observation, and each column represents
#'             one predictor
#'
#'@param  lambda: positive penalty parameter value
#'@param      b0: numeric initial start of the regression vector
#'@param    sig0: numeric positive scalar, initial scale estimate.
#'@param  c: Threshold constant of Huber's loss function
#'@param  reltol: Convergence threshold. Terminate when successive
#'             estimates differ in L2 norm by a rel. amount less than reltol.
#'             Default is 1.0e-5
#'@param  iter_max: int, default = 500. maximum number of iterations
#'@param  printitn: print iteration number (default = 0, no printing)
#'
#'@return        b0: regression coefficient vector estimate
#'@return     sig0: estimate of the scale
#'@return   psires: pseudoresiduals
#'
#'@examples
#'hublasso(rnorm(5), matrix(rnorm(5)), lambda = 0.5, b0 = rnorm(5), sig0 = 0.3)
#'
#'@note
#'
#'File in Regression.R
#'
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
  normb0 <- norm(b0, type = "2")

  for(iter in 1:iter_max){
    r <- y - X %*% b0
    psires <- psihub(r/sig0, c) * sig0
    # scale update
    sig1 <- norm(psires, type = "2") / con

    crit2 <- abs(sig0 - sig1)
    for(jj in 1:p){
      # update pseudoresidual
      psires <- psihub(r/sig1, c) * sig1
      b0[jj] <- soft_thresh(b0[jj] + X[,jj] %*% psires, lambda)
      r <- r + X[,jj] * (betaold[jj] - b0[jj])
    }

    normb <- norm(b0, type = "2")

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
            lambda,iter, norm(fpeq, type = "2"),sig1-norm(psires, type = "2")/con)
  }

  return (list('b0' = b0,'sig0' = sig0, 'psires' = psires))
}






# ladlassopath ----

#' ladlassopath
#'
#' ladlassopath computes the LAD-Lasso regularization path (over grid
#'                                                           of penalty parameter values). Uses IRWLS algorithm.
#'
#' @param       y : Numeric data vector of size N x 1 (output, respones)
#' @param     X : Numeric data matrix of size N x p. Each row represents one
#'          observation, and each column represents one predictor (feature).
#' @param intcpt: Logical (true/false) flag to indicate if intercept is in the
#'          regression model
#' @param eps: Positive scalar, the ratio of the smallest to the
#'          largest Lambda value in the grid. Default is eps = 10^-3.
#' @param     L : Positive integer, the number of lambda values EN/Lasso uses.
#'          Default is L=120.
#' @param reltol : Convergence threshold for IRWLS. Terminate when successive
#'          estimates differ in L2 norm by a rel. amount less than reltol.
#' @param printitn: print iteration number (default = 0, no printing)
#'
#' @return   B    : Fitted LAD-Lasso regression coefficients, a p-by-(L+1) matrix,
#'         where p is the number of predictors (columns) in X, and L is
#'         the  number of Lambda values. If intercept is in the model, then
#'         B is (p+1)-by-(L+1) matrix, with first element the intercept.
#' @return stats  : structure with following fields:
#'            Lambda = lambda parameters in ascending order
#'            MeAD = Mean Absolute Deviation (MeAD) of the residuals
#'            gBIC = generalized Bayesian information criterion (gBIC) value
#'                  for each lambda parameter on the grid.
#'@examples
#'ladlassopath(rnorm(5), matrix(rnorm(5)))
#'
#' @note
#' File in Regression.R
#'
#' @export
ladlassopath <- function(y, X, L = 120, eps = 1e-3, intcpt = T, reltol = 1e-6, printitn = 0){
  n = nrow(X)
  p = ncol(X)

  if(intcpt){
    p <- p + 1
    if(is.complex(y)) medy <- spatmed(y) else medy <- median(y) # spatmed is from 03_covariance
    yc <- y - medy
    lam0 <- norm(t(X) %*% sign(yc), type = "I")
    } else lam0 <- norm(t(X) %*% sign(y), type = "I")

  lamgrid <- eps^(0:L / L) * lam0
  B <- diag(0, p, L + 1)

  if(intcpt) binit <- c(medy, numeric(p-1)) else binit <- numeric(p)

  for(jj in 1:(L+1)){
    B[,jj] <- ladlasso(y,X, lamgrid[jj], intcpt, binit, reltol, printitn)[[1]]
    binit <- B[,jj]
  }

  if(intcpt){
    B[rbind(rep(F, L+1), abs(B[2:nrow(B),]) < 1e-7)] <- 0
    DF <- colSums(abs(B[2:nrow(B),,drop = F]) != 0)
    MeAD <- sqrt(pi / 2) * colMeans(abs(repmat(y, L+1) - cbind(rep(1, n), X) %*% B))
    const <- sqrt(n / (n - DF - 1))
  } else {
    B[abs(B) < 1e-7] <- 0
    DF <- colSums(abs(B) != 0)
    MeAD <- sqrt(pi / 2) * colMeans(abs(repmat(y, L+1) - X %*% B))
    const <- sqrt(n / (n - DF - 1))
  }

  gBIC <- 2 * n * log(MeAD) + DF * log(n)

  stats <- list(MeAD * const, gBIC, lamgrid, names = c('MeAD', 'gBIC', 'Lambda'))

  return(list(B, stats))
}

# ladreg ----

#' ladreg
#'
#' ladreg computes the LAD regression estimate
#'
#' @param          y: numeric response N vector (real/complex)
#' @param        X: numeric feature  N x p matrix (real/complex)
#' @param   intcpt: (logical) flag to indicate if intercept is in the model
#' @param       b0: numeric optional initial start of the regression vector for
#' IRWLS algorithm. If not given, we use LSE (when p>1).
#' @param printitn: print iteration number (default = 0, no printing) and
#' other details
#'
#' @return        b1: (numberic) the regression coefficient vector
#' @return     iter: (numeric) # of iterations (given when IRWLS algorithm is used)
#'
#' @examples
#'ladreg(rnorm(5), matrix(rnorm(5)))
#'@note
#'file location: Regression.R
#' @export
ladreg <- function(y, X, intcpt = T, b0 = NULL, printitn = 0){
  return(ladlasso(y, X, 0, intcpt, b0, printitn))

}

# Mreg ----

#' Mreg computes the M-estimates of regression using an auxiliary scale
#' estimate. It uses the iterative reweighted least squares (IRWLS) algorithm
#'
#' @param        y : (numeric) data vector of size N (output, response vector)
#' @param      X : (numeric) data matrix of size N x p (input, feature matrix)
#'           If the model has intercept, then first column of X should be a
#'           vector of ones.
#' @param lossfun : (string) either 'huber' or 'tukey' to identify the desired
#'           loss function. Default is 'huber'
#' @param     b0 : (numeric) Optional robust initial start (regression vector) of
#'           iterations. If not given, we use the LAD regression estimate
#' @param verbose: (logical) true of false (default). Set as true if you wish
#'           to see convergence as iterations evolve.
#'
#' @return b1 : regression parameters
#' @return sig: scale
#'
#' @examples
#'Mreg(1:5, matrix(-1:3))
#'Mreg(1:5, matrix(-1:3), lossfun = 'tukey')
#' @export
Mreg <- function(y, X, lossfun = 'huber', b0 = NULL, verbose = F){

  if(is.null(b0)) b0 <- ladreg(y, X, F)$b1[[1]]

  # Compute the auxiliary scale estimate as
  if(is.complex(y)) const <- 1.20112 else const <- 1.4815

  resid <- abs(y - X %*% b0) # resid is a matrix

  sig <- const * median(resid[resid !=  0]) # auxiliary scale estimate

  if(lossfun == 'tukey'){
    if(is.complex(y)) c <- 3 else c <- 3.4437
    wfun <- function(x) wtuk(x, c)
  }
  else if(lossfun == 'huber'){
    if(is.complex(y)) c <- 1.214 else c <- 1.345
    wfun <- function(x) whub(x, c)
  } else {stop('lossfun must be either huber or tukey as a string')}

  ITERMAX <- 1000
  TOL <- 1.0e-5

  if(verbose) sprintf('Mreg: iterations starting, using %s loss function \n',lossfun)

  for(iter in 1:ITERMAX){

    resid[resid < 1e-6] <- 1e-6
    w <- wfun(resid / sig)
    Xstar <- X * w
    b1 <- (t(Xstar) %*% y) / (t(Xstar) %*% X)

    crit <- norm(b1 - b0, type = "2") / norm(b0, type = "2")
    if(verbose) sprintf('Mreg: crit(%4d) = %.9f\n',iter,crit)
    if(crit < TOL) break

    b0 <- b1

    resid <- abs(y - X %*% b0)
  }

  return(list('b1' = b1, 'sig' = sig))
  }

# rankflasso ----

#' Computes the rank fused-Lasso regression estimates for given fused
#' penalty value lambda_2 and for a range of lambda_1 values
#'
#' @param  y       : numeric response N x 1 vector (real/complex)
#' @param  X       : numeric feature  N x p matrix (real/complex)
#' @param  lambda1 : positive penalty parameter for the Lasso penalty term
#' @param  lambda2 : positive penalty parameter for the fused Lasso penalty term
#' @param  b0      : numeric optional initial start (regression vector) of
#'                   iterations. If not given, we use LSE (when p>1).
#' @param printitn : print iteration number (default = 0, no printing)
#'
#' @return    b      : numeric regression coefficient vector
#' @return  iter   : positive integer, the number of iterations of IRWLS algorithm
#'
#' @examples
#'
#' @export
rankflasso <- function(y, X, lambda1, lambda2, b0 = NULL, printitn = 0){
  n <- nrow(X)
  p <- ncol(X)

  intcpt <- F

  if(is.null(b0)) {
    b0 <- qr.solve(cbind(rep(1, n), X), y)
    b0 <- b0[2:length(b0)]
  }

  B <- repmat(1:n, n)

  A <- t(B)

  a <- A[A < B]
  b <- B[A < B]

  D <- sparseMatrix(1:(p-1), 1:(p-1), x = -1, dims = c(399, 400)) # diag(-1, p-1, p-1)
  D[seq(p, p*(p - 1), p)] <- 1

  ytilde <- sparseVector(x = y[a] - y[b], i = 1:length(a), length = length(a)+ p - 1)

  Xtilde <- rbind( X[a,] - X[b,], lambda2 * D) # 80199 400

  if(printitn > 0) sprintf('rankflasso: starting iterations\n')

  r <- ladlasso(ytilde, Xtilde, lambda1, intcpt, b0, printitn)
  iter <- r[[2]]
  b <- matrix(r[[1]]) # wrapped in matrix
  b[abs(b) < 1e-7] <- 0
  return(list(b,iter))
}

# rankflassopath ----

#' rankflassopath()
#'
#' Computes the rank fused-Lasso regression estimates for given fused
#' penalty value lambda_2 and for a range of lambda_1 values
#'
#' @param    y       : numeric response N x 1 vector (real/complex)
#' @param    X       : numeric feature  N x p matrix (real/complex)
#' @param   lambda2 : positive penalty parameter for the fused Lasso penalty term
#' @param   L       : number of grid points for lambda1 (Lasso penalty)
#' @param   eps     : Positive scalar, the ratio of the smallest to the
#'             largest Lambda value in the grid. Default is eps = 10^-4.
#' @param printitn : print iteration number (default = F, no printing)
#'
#'
#'@return         B: Fitted rank fused-Lasso regression coefficients, a p-by-(L+1) matrix,
#'           where p is the number of predictors (columns) in X, and L is
#'           the  number of Lambda values.
#'@return       B0: estimates values of intercepts
#'@return       lamgrid: = lambda parameters
#'@examples
#'
#'@export
rankflassopath <- function(y, X, lambda2, L = 120, eps = 1e-3, printitn = F){
  n <- nrow(X)
  p <- ncol(X)

  intcpt <- F

  B <- repmat(1:n, n)
  A <- t(B)

  a <- A[A < B]
  b <- B[A < B]

  D <- diag(-1, p-1, p-1)
  D[seq(p, (p-1)^2, p)] <- 1

  onev <- c(rep(0, p-2), 1)
  D <- cbind(D, onev)

  ytilde <- c(y[a] - y[b], rep(0, p-1))

  Xtilde <- rbind( X[a,] - X[b,], lambda2 * D)

  lam0 <- norm(t(Xtilde) %*% mat_sign(ytilde), type = "I")

  lamgrid <- eps^(0:L / L) * lam0

  B <- rep(0, L + 1)
  B0 <- rep(0, L + 1)

  b_init <- rep(0, p)

  if(printitn) sprintf('rankflassopath: starting iterations\n')

  for(jj in 1:(L+1)){
    B[, jj] <- ladlasso(ytilde, Xtilde, lamgrid[jj], intcpt, b_init, printitn)[[1]]
    b_init <- B[, jj]
    if(printitn) print('.')
  }

  B[abs(B) < 1e-7] = 0

  result <- list(B, B0, lamgrid)

  names(result) <- c('B', 'B0', 'lamgrid')

  return( result)
}

# ranklasso ----

#' ranklasso
#'
#' ranklasso computes the rank (LAD-)regression estimate
#'
#' @param        y  : numeric data vector of size N x 1 (output, respones)
#' @param      X  : numeric data matrix of size N x p (input, features)
#' @param  lambda : penalty parameter (>= 0)
#' @param     b0  : numeric optional initial start (regression vector) of
#'           iterations. If not given, we use LSE.
#' @param printitn : bool, print iteration number (default = F, no printing) and
#'            other details
#'
#' @return    b1     : numeric the regression coefficient vector
#' @return   iter   : (numeric) # of iterations (given when IRWLS algorithm is used)
#'
#' @examples
#'
#' @export
ranklasso <- function(y, X, lambda, b0 = NULL, printitn = F){
  n <- nrow(X)

  intcpt <- F

  if(is.null(b0)){
    b0 <- qr.solve(cbind(rep(1, n), X), y)
    b0 <- b0[2:length(b0)]
  }

  B <- repmat(1:n, n)

  A <- t(B)

  a <- A[A < B]
  b <- B[A < B]

  Xtilde <- X[a,] - X[b,]
  ytilde <- y[a] - y[b]

  return(ladlasso(ytilde, Xtilde, lambda, intcpt, b0, printitn))
}


# rladreg ----

#' rladreg
#'
#'computes the LAD regression estimate
#'
#' @param        y: numeric response N vector (real/complex)
#' @param        X: numeric feature  N x p matrix (real/complex)
#' @param       b0: numeric optional initial start of the regression vector for
#'                  IRWLS algorithm. If not given, we use LSE (when p>1).
#' @param printitn: print iteration number (default = 0, no printing) and
#'            other details
#'
#' @return          b1: numeric the regression coefficient vector
#' @return        iter: (numeric) # of iterations (given when IRWLS algorithm is used)
#'
#'
#' @examples
#' library(MASS)
#' rladreg(1:5, matrix(-1:3))
#' rladreg(1:5 +1i, matrix(-1:3))
#' @note
#' file location: Regression.R
#' @export
rladreg <- function(y, X, b0 = NULL, printitn = F){
  n <- nrow(X)
  p <- ncol(X)

  if(is.null(b0)){
    b0 <- ginv(cbind(rep(1, n), X)) %*% y
    b0 <- b0[2:length(b0)]
  }

  B <- repmat(1:n, n)
  A <- t(B)
  a <- A[A < B]
  b <- B[A < B]
  Xtilde <- X[a,] - X[b,]
  ytilde <- y[a] - y[b]

  if(p == 1){
    b1 <- wmed(ytilde / Xtilde, abs(Xtilde))[[1]]
    iter <- NULL
    return(list('b1' = b1, 'iter' = iter))
  } else return(ladlasso(ytilde, Xtilde, 0, b0, printitn))
}

# wmed ----

#' weighted median
#'
#' wmed computes the weighted median for data vector y and weights w, i.e.
#' it solves the optimization problem:
#'
#' beta = arg min_b  SUM_i | y_i - b | * w_i
#'
#'
#'@param        y : (numeric) data given (real or complex)
#'@param       w : (nubmer) positive real-valued weights. Inputs need to be of
#'           same length
#'@param       verbose: (logical) true of false (default). Set as true if you wish
#'       to see convergence as iterations evolve
#'
#'@param tol: threshold for iteration in complex case \cr default = 1e-7
#'@param iter_max: number of iterations in complex case \cr default = 2000
#'
#'@return beta: (numeric) weighted median
#'@return converged: (logical) flag of convergence (in the complex-valued
#'                                                         data case)
#'@return iter: (numeric) the number of iterations (complex-valued case)
#'
#'@examples
#'wmed(1:5, c(1,0,1,2,3))
#'wmed(1:5 +1i, c(1,0,1,2,3))
#'
#'@export
wmed <- function(y, w, verbose = F, tol = 1e-7, iter_max = 2000){
  N <- length(y)

  # complex-valued case
  if(is.complex(y)){
    beta0 <- median(Re(y)) + 1i * median(Im(y)) # initial guess
    abs0 <- abs(beta0)
    for(iter in 1:iter_max){
      wy <- abs(y - beta0)
      wy[wy <= 1e-6] <- 1e-6
      update <- sum(w * mat_sign(y - beta0)) / sum(w / wy)
      beta <- beta0 + update
      delta <- abs(update) / abs0
      if(verbose & iter %% 10 == 0) sprintf('At iter = %3d, delta=%.8f\n',iter,delta)
      if(delta <= tol) break

      beta0 <- beta
      abs0 <- abs(beta)
    }# ends for loop
    if(iter == iter_max) converged <- F else converged <- T
  } else{
    # real- value case
    tmp <- sort(y, index.return = T)
    y <- tmp[[1]]
    indx <- tmp[[2]]

    w <- w[indx]
    wcum <- cumsum(w)
    i <- which(wcum < 0.5 * sum(w))
    i <- i[length(i)]
    beta <- y[i+1] # see equation (2.21) of the book

    iter <- 0
    converged <- NULL
  } # end real value case

  return( list('beta' = beta, 'iter' = iter, 'converged' = converged))
}
