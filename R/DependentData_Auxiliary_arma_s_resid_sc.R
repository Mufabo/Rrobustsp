arma_s_resid_sc <- function(x, beta_hat, p, q) {
  # phi_hat := AR coefficients
  if (0 < p) { phi_hat <- beta_hat[1:p]
  } else phi_hat <- numeric(0)
  if (0 < q) { theta_hat <- beta_hat[(p + 1):(p + q)]
  } else theta_hat <- numeric(0)

  N <- length(x)
  r <- max(p, q)

  a <- numeric(N)
  x_sc <- m_scale(x)

  if (r == 0) {
    a <- x
    a_sc <- x_sc
  }
  else{
    if (or(sum(abs(roots(c(1, -phi_hat))) > 1), sum(abs(roots(c(1, -theta_hat))) > 1))){
      return(10 ^ 10)
      }
    else{
      #ARMA
      if (p >= 1 & q >= 1) {
        for (ii in (r + 1):N) {
          # ARMA residuals
          a[ii] <-
            x[ii] - phi_hat %*% x[seq.int(ii - 1, ii - p, by = -1)] + theta_hat %*% a[seq.int(ii -
                                                                                                1, ii - q, by = -1)]
        }
      }
      # MA model
      if (p == 0 & q >= 1) {
        for (ii in (r + 1):N) {
          a[ii] <- x[ii] + theta_hat %*% a[seq.int(ii - 1, ii - q, by = -1)]
        }
      }
      # AR model
      if (p >= 1 & q == 0) {
        for (ii in (r + 1):N) {
          a[ii] <- x[ii] - phi_hat %*% x[seq.int(ii - 1, ii - p, by = -1)]
        }
      }

    }
  }

  a_sc <- m_scale(a[(p + 1):N])
  return(a_sc)
}
