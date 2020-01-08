#' bip_ar1_s
#'
#' @export
bip_ar1_s <- function(x, N, phi_grid, fine_grid, kap2) {
  # AR(0): residual scale equals observation scale
  a_scale_final <- m_scale(x) # ok

  # grid search for partial autocorrelations
  a_bip_sc <- numeric(length(phi_grid))
  a_sc <- numeric(length(phi_grid))

  for (mm in 1:length(phi_grid)) {
    a <- numeric(length(x)) # residuals for BIP-AR

    a2 <- numeric(length(x)) # residuals for AR

    lambda <- ma_infinity(phi_grid[mm], 0, 100) # ok

    # sigma used for BIP-model
    sigma_hat <- a_scale_final / sqrt(1 + kap2 * sum(lambda ^ 2)) # ok

    for (ii in 2:N) {
      # residuals for BIP-AR
      a[ii] <- x[ii] - phi_grid[mm] * (
        x[ii - 1] - a[ii - 1] + sigma_hat * Rrobustsp::eta(a[ii - 1] / sigma_hat))

      # residuals for AR
      a2[ii] <- x[ii] - phi_grid[mm] * x[ii - 1]
    }

    # tau-scale of residuals for BIP-AR
    a_bip_sc[mm] <- m_scale(a[2:N])

    # tau-scale of residuals for AR
    a_sc[mm] <- m_scale(a2[2:N])
  } # for mm in 1:length(phi_grid)

  # polynomial approximation of tau scale objective function for BIP-AR(1) tau-estimates
  poly_approx <- pracma::polyfit(phi_grid, a_bip_sc, 5)

  # interpolation of tau scale objective function for BIP-AR(1) tau-estimates to fine grid
  a_interp_scale <- polyval(poly_approx, fine_grid)

  # polynomial approximation of  tau scale objective function for AR(1) tau-estimates
  poly_approx2 <- polyfit(phi_grid, a_sc, 5)

  # interpolation of  tau scale objective function for AR(1) tau-estimates to fine grid
  a_interp_scale2 <- polyval(poly_approx2, fine_grid)

  temp <- min(a_interp_scale)
  ind_max <- which.min(a_interp_scale)

  # tau-estimate under the BIP-AR(1)
  phi <- fine_grid[ind_max]

  temp2 <- min(a_interp_scale2)
  ind_max2 <- which.min(a_interp_scale2)

  # tau-estimate under the AR(1)
  phi2 <- fine_grid[ind_max2]

  # final estimate maximizes robust likelihood of the two
  if (temp2 < temp) {
    phi_s <- phi2
    temp <- temp2
  } else {
    phi_s <- phi
  }

  # final BIP-tau-estimate for AR(1)
  phi_hat <- phi_s

  # final AR(1) tau-scale-estimate depending on phi_hat
  lambda <- ma_infinity(phi_hat, 0, 100)

  sigma_hat <- a_scale_final / sqrt(1 + kap2 * sum(lambda ^ 2))

  a <- numeric(length(x))
  a2 <- numeric(length(x))
  x_filt <- numeric(length(x))

  for (ii in 2:N) {
    a[ii] <- x[ii] - phi_hat * (x[ii - 1] - a[ii - 1] + sigma_hat *
                                  Rrobustsp::eta(a[ii - 1] / sigma_hat))
    x_filt[ii] <- x[ii] - a[ii] + sigma_hat * Rrobustsp::eta(a[ii] / sigma_hat)

    a2[ii] <- x[ii] - phi_hat * x[ii - 1]
  }

  if(temp2 < temp) a_scale_final <- c(a_scale_final, m_scale(a[2:N])) else a_scale_final <- c(a_scale_final, m_scale(a2[2:N]))

  return(list(
    'x_filt' = x_filt,
    'phi_hat' = phi_hat,
    'a_scale_final' = a_scale_final
  ))
}
