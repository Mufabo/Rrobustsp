#' ar_est_bip_tau
#'
#' The function ar_est_bip_tau(x,P) comuptes BIP S-estimates of the
#' AR model parameters. It also computes an outlier cleaned signal using
#' BIP-AR(P) predictions, and the M-scale of the estimated innovations
#' series.
#'
#' @param x: data (observations/measurements/signal)
#' @param P: autoregressive order
#'
#' @return phi_hat: a vector of bip tau estimates for each order up to P
#' @return x_filt: cleaned version of x using robust BIP predictions
#' @return a_scale_final: minimal tau scale of the innovations of BIP AR or AR
#'
#' @examples
#'
#' @references
#'
#' "Robust Statistics for Signal Processing"
#' Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#' Cambridge University Press, 2018.
#'
#'"Bounded Influence Propagation \eqn{\tau}-Estimation: A New Robust Method for ARMA Model Estimation."
#' Muma, M. and Zoubir, A.M.
#' IEEE Transactions on Signal Processing, 65(7), 1712-1727, 2017.
#' @export
ar_est_bip_tau <- function(x, P) {
  N <- length(x)
  kap2 = 0.8724286 # kap=var(eta(randn(10000,1)))

  # coarse grid search
  phi_grid <- seq(-.99, .99, by = 0.05)

  # finer grid via polynomial interpolation
  fine_grid <- seq(-.99, .99, by = 0.001)

  a_bip_sc <- numeric(length(phi_grid))
  a_sc <- numeric(length(phi_grid))

  # The following was introduced so as not to predict based on highly contaminated data in the
  # first few samples.
  x_tran <- x[1:min(10, floor(N / 2))]
  sig_x_tran <- madn(x)
  x_tran[abs(x_tran) > 3 * sig_x_tran] <-
    sign(abs(x_tran) > 3 * sig_x_tran) * 3 * sig_x_tran
  x[1:min(10, floor(N / 2))] <- x_tran[1:min(10, floor(N / 2))]

  a_scale_final <- numeric(P + 1)

  # P == 0 ----
  if (P == 0) {
    phi_hat <- NULL
    # AR(0) residual scale equals observation scale
    a_scale_final[1] <- tau_scale(x)

  }

  # P == 1 ----
  if (P == 1) {
    # AR(1) tau-estimates
    c(x_filt, phi_hat, a_scale_final) %<-% bip_ar1_tau(x, N, phi_grid, fine_grid, kap2)
  }

  # P > 1 ----
  if (P > 1) {
    phi_hat <- matrix(0, P, P)

    # AR(1) tau-estimates
    c(x_filt, phi_hat[1, 1], a_scale_final) %<-% bip_ar1_tau(x, N, phi_grid, fine_grid, kap2)

    for (p in 2:P) {
      for (mm in 1:length(phi_grid)) {
        for (pp in 1:(p - 1)) {
          phi_hat[p, pp] <-
            phi_hat[p - 1, pp] - phi_grid[mm] * phi_hat[p - 1, p - pp]
        }

        predictor_coeffs <- c(phi_hat[p, 1:(p - 1)], phi_grid[mm])
        M <- length(predictor_coeffs)

        if (mean(abs(roots(c(
          1, predictor_coeffs
        ))) < 1) == 1) {
          lambda <- ma_infinity(predictor_coeffs, 0, 100)
          # sigma used for bip-model
          sigma_hat <-
            a_scale_final[1] / sqrt(1 + kap2 * sum(lambda ^ 2))
        } else {
          sigma_hat <- mad(x, constant = 1.483)
        }

        a <- numeric(length(x))
        a2 <- numeric(length(x))

        for (ii in (p + 1):N) {
          intv <- seq.int(ii - 1, ii - M, by = -1)
          a[ii] <- x[ii] - predictor_coeffs %*%
            (x[intv] - a[intv] + sigma_hat * Rrobustsp::eta(a[intv] / sigma_hat))

          a2[ii] <- x[ii] - predictor_coeffs %*% x[intv]
        }

        # residual scale for BIP-AR
        a_bip_sc[mm] <- tau_scale(a[(p + 1):length(a)])
        # residual scale for AR
        a_sc[mm] <- tau_scale(a2[(p + 1):length(a2)])
      } # for (mm in 1:length(phi_grid))

      c(ind_max, temp, ind_max2, temp2) %<-% res_scale_approx(phi_grid, a_bip_sc, fine_grid, a_sc)

      # tau-estimate under the BIP-AR(p)
      phi <- -fine_grid[ind_max]

      # tau-estimate under the AR(p)
      phi2 <- -fine_grid[ind_max2]

      # final estimate minimizes the residual scale of the two
      if (temp2 < temp) {
        ind_max <- ind_max2
        temp <- temp2
      }

      for (pp in 1:(p - 1)) {
        phi_hat[p, pp] <- phi_hat[p - 1, pp] - fine_grid[ind_max] *
          phi_hat[p - 1, p - pp]
      }

      phi_hat[p, p] <- fine_grid[ind_max]

      # final AR(P) tau-scale-estimate depending on phi_hat(p,p)
      if (mean(abs(roots(c(1, phi_hat[p,]))) < 1) == 1) {
        lambda <- ma_infinity(phi_hat[p, ], 0, 100)
        # sigma used for bip-model
        sigma_hat <- a_scale_final[1] /
          sqrt(1 + kap2 * sum(lambda ^ 2))
      } else
        sigma_hat <- mad(x, constant = 1.483)
      # hier phi-hat wrong
      x_filt <- numeric(length(x))

      for (ii in (p + 1):N) {
        intv <- seq.int(ii - 1, ii - M, by = -1)
        a[ii] <- x[ii] - phi_hat[p, 1:p] %*%
          (x[intv] - a[intv] +
             sigma_hat * Rrobustsp::eta(a[intv] / sigma_hat))

        a2[ii] <- x[ii] - phi_hat[p, 1:p] %*% x[intv]
      }

      if (temp < temp2)
        a_scale_final[p + 1] <- tau_scale(a[(p + 1):N])
      else
        a_scale_final[p + 1] <- tau_scale(a2[(p + 1):N])
    } # for p in 2:P

    # BIP-AR(P) tau-estimates
    phi_hat <- phi_hat[p, ]

    for (ii in (p + 1):N) {
      x_filt[ii] <-
        x[ii] - a[ii] + sigma_hat * Rrobustsp::eta(a[ii] / sigma_hat)
    }
  }

  return(list(
    'phi_hat' = phi_hat,
    'x_filt' = x_filt,
    'a_scale_final' = a_scale_final
  ))
}
