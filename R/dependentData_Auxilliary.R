
#' eta
#'
#' @param x: numeric vector, the signal
#' @param c: numeric, default = 1
#'
#' @examples
#' library(Rrobustsp)
#'
#' x <- rnorm(5)
#' eta(x)
#'
#' @export
eta <- function(x, c = 1){
  x <- x / c

  y <- x

  y[abs(x) > 3] <- 0

  y[abs(x) > 2 & abs(x) <= 3] <- 0.016 * x[abs(x) > 2 & abs(x) <= 3]^7 -
   0.312 * x[abs(x) > 2 & abs(x) <= 3]^5 +
   1.728 * x[abs(x) > 2 & abs(x) <= 3]^3 -
   1.944 * x[abs(x) > 2 & abs(x) <= 3]

  y[abs(x) <= 2] <- x[abs(x) <= 2]

  y <- c * y
  return(y)
}

ma_infinity <- function(phi, theta, Q_long){
  Q <- length(theta)
  P <- length(phi)
  theta_inf <- pracma::deconv(c(1, theta, numeric(Q_long + P + Q)), c(1, -phi))$q
  theta_inf <- theta_inf[2:(Q_long + 1)]
  return(theta_inf)
}

muler_rho1 <- function(x){
  x <- x / 0.405 # where does the 0.405 come from ?

  intv <- abs(x) > 2 & abs(x) <= 3

  rho <- numeric(length(x))

  rho[abs(x) <= 2] <- 0.5 * x[abs(x) <= 2]^2
  rho[intv] <- 0.002 * x[intv]^8 - 0.052 * x[intv]^6 + 0.432 * x[intv]^4 - 0.972 * x[intv]^2 + 1.792

  rho[abs(x) > 3] <- 3.25
  return(rho)
}


#' muler_rho2
#'
#'
#' @note
#' Location: .../Rrobustsp/dependentData_Auxiliary.R
#'
#' @export
muler_rho2 <- function(x){
  rho <- rep(3.25, length(x))

  intv <- (abs(x) > 2) & (abs(x) <= 3)
  rho[intv] <-0.002 * x[intv]^8- 0.052 * x[intv]^6+ 0.432 * x[intv]^4- 0.972 * x[intv]^2+ 1.792

  rho[abs(x) <= 2] <- 0.5 * x[abs(x) <= 2]^2
  return(rho)
}

m_scale <- function(x){
  N <- length(x)
  sigma_k <- madn(x)
  delta <- 3.25 / 2 # max(muler_rho1)/2

  epsilon <- 1e-4
  w_k <- rep(1, N)

  max_iters <- 30
  k <- 0

  while(k<=max_iters & sigma_k < 10^5){
    w_k[x != 0] <- muler_rho1(x[x != 0] / sigma_k) / (x[x != 0] / sigma_k)^2

    sigma_k_plus1 <- sqrt(1 / (N * delta) * sum(w_k * x^2))
    if(!is.nan(sigma_k_plus1 / sigma_k -1) & abs(sigma_k_plus1 / sigma_k -1) > epsilon){
      sigma_k <- sigma_k_plus1
      k <- k + 1
    } else break
  }
  sigma_hat <- sigma_k
  return(sigma_hat)
}

res_scale_approx <- function(phi_grid, a_bip_sc, fine_grid, a_sc) {
  # polynomial approximation of residual scale for BIP-AR(p) tau-estimates
  poly_approx <- polyfit(phi_grid, a_bip_sc, 5)

  # interpolation of  residual scale for BIP-AR(p) tau-estimates to fine grid
  a_interp_scale <- c(polyval(poly_approx, fine_grid))

  # polynomial approximation of  residual scale for AR(p) tau-estimates
  poly_approx2 <- polyfit(phi_grid, a_sc, 5)

  # interpolation of  residual scale for AR(p) tau-estimates to fine grid
  a_interp_scale2 <- c(polyval(poly_approx2, fine_grid))

  temp <- min(a_interp_scale)
  ind_max <- which.min(a_interp_scale)

  temp2 <- min(a_interp_scale2)
  ind_max2 <- which.min(a_interp_scale2)

  return(list('ind1' = ind_max, 'min1' = temp, 'ind2' = ind_max2, 'min2' = temp2))
}

#' tau_scale
#'
#'
#' @param x
#'
#' @return scale
#'
#' @note
#' Location: .../Rrobustsp/R/dependentData_Auxiliary
#'
#' @export
tau_scale <- function(x){
  b <- 0.398545548533895;  # E(muler_rho2) under the standard normal distribution
  sigma_m <- m_scale(x);
  sigma_hat <- sqrt(sigma_m^2/(length(x)) * 1/b * sum(muler_rho2(x/sigma_m)));
  return(sigma_hat)
}
