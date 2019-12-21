# Args ----

X <- diag(1, 400, 400)
load('~/Rrobustsp/data/images.RData')
load(path_test('blas'))

y <- unlist(unname(images['y20n']))
lambda2 <- 340
lambda1 <- 124
b0 <- Blas

printitn <- 1

# ----
n <- nrow(X)
p <- ncol(X)

intcpt <- F

if(is.null(b0)) {
  b0 <- MASS::ginv(cbind(rep(1, n), X)) %*% y #qr.solve(cbind(rep(1, n), X), y)
  b0 <- bo[2:length(b0)]
}

B <- repmat(1:n, n)

A <- t(B)

a <- A[A < B]
b <- B[A < B]

D <- sparseMatrix(1:(p-1), 1:(p-1), x = -1) # diag(-1, p-1, p-1)
D[2:(p-1), 2:(p-1)] <- 1
D <- cbind(D, c(rep(0, p-2), 1)) # 399 400

ytilde <- sparseVector(x = y[a] - y[b], i = 1:length(a), length = length(a)+ p - 1)

Xtilde <- rbind( X[a,] - X[b,], lambda2 * D) # 80199 400

if(printitn > 0) sprintf('rankflasso: starting iterations\n')

# args ladlasso ----

# r <- ladlasso(ytilde, Xtilde, lambda1, intcpt, b0, reltol = 1)
y <- ytilde
X <- Xtilde
lambda <- lambda1
iter_max <- 100
reltol <- 1
# ladlasso ----

  N <- nrow(X)
  p <- ncol(X)

  if(intcpt) X <- cbind(matrix(1, N, 1), X)

  if(is.null(b0)) b0 <- qr.solve(X, y) # ginv(X) %*% y

  iter <- NULL
  if(printitn > 0) sprintf('Computing the solution for lambda = %.3f\n',lambda)

  if(printitn > 0) print('Starting the IRWLS algorithm..\n')

  if(lambda > 0){
    y <- c(y, rep(0, p))

    # slow
    if(intcpt) X <- rbind(X, cbind(rep(0, p), diag(lambda, p, p))) else X <- rbind(X, diag(lambda, p, p))
  }

  for(iter in 1:iter_max){
    resid <- abs(y - X %*% b0)
    resid[resid < 1e-6] <- 1e-6

    # slow
    Xstar <- sweep(X, 1, resid, FUN = "/")

    # slow
    # b1 <- fRcpp(Xstar, X, matrix(y))

    b1 <- qr.solve(t(Xstar) %*% X, t(Xstar) %*% y)

    # b1 <- c(MASS::ginv(t(Xstar) %*% X) %*% (t(Xstar) %*% y))

    crit <- norm(b1-b0, type = "2") / norm(b0, type = "2")

    if(printitn > 0 & iter %% printitn) sprintf('ladlasso: crit(%4d) = %.9f\n',iter,crit)
    if(crit < reltol && iter > 10) break
    b0 <- b1
  }

# ladlasso end ----

b1[abs(b1) < 1e-7] <- 0
expect_equal(iter, 11)

load(path_test('reg_rankflasso_1'))
expect_equal(round(b1, digits = 3), sol)
