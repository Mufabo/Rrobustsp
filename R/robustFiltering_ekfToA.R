#' ekf_toa
#'
#' EKF for tracking with time-of-arrival (:= ToA) estimates
#'
#' @param r_ges : measured distances as a M x N matrix
#' @param theta_init : initial state estimate
#' @param BS : base station positions
#' @param parameter :
#'
#' @return th_hat : state estimates
#' @return P_min : apriori covariance
#' @return P : aposteriors covariance
#'
#' @examples
#' ekf_toa(matrix(rnorm(16), 4, 4),
#'         rnorm(4),
#'         matrix(rnorm(8), 4, 2))
#' @references
#' "Robust Statistics for Signal Processing"
#'  Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#'  Cambridge University Press, 2018.
#'
#'  "Robust Tracking and Geolocation for Wireless Networks in NLOS Environments."
#'   Hammes, U., Wolsztynski, E., and Zoubir, A.M.
#'   IEEE Journal on Selected Topics in Signal Processing, 3(5), 889-901, 2009.
#' @export
ekf_toa <- function(r_ges, theta_init, BS, parameter = NULL){
  if(is.null(parameter)){
    print('Using default parameters')
    sigma_v <- 1
    M <- nrow(BS)
    P0 <- diag(x = c(100, 100, 10, 10))
    R <- 150^2 * diag(x = 1, M, M)
    Ts <- 0.2
    A <- matrix(data = c(1.,0,0,0, 0,1,0,0, Ts,0,1,0, 0,Ts,0,1), 4, 4)
    Q <- sigma_v^2 * diag(1,2,2)
    G <- rbind(Ts^2 / 2 * diag(1, 2, 2), Ts * diag(1,2,2))
    #dimension of positions, default is 2
    parameter <- list('dim' = ncol(BS), 'var.est' = 1)
  } else {
    P0 <- parameter$P0
    R <- parameter$R
    Q <- parameter$Q
    G <- parameter$G
    A <- parameter$A
  }

  if(2 * parameter$dim != length(theta_init) |
     2 * parameter$dim != nrow(P0)){
    simpleError('State vector or state covariance do not match the dimensions of the BS')
  }

  x <- BS[, 1]
  y <- BS[, 2]

  M <- length(x)
  N <- ncol(r_ges)

  P <- to.tensor(0,c(U=4,V=4,W=N))
  th_hat <- matrix(theta_init, length(theta_init), N)
  th_hat_min <- matrix(0, 4, N)
  P_min <- to.tensor(0,c(U=4,V=4,W=N))
  H <- matrix(0, M, 4)
  h_min <- numeric(M)
  sigma2 <- numeric(N)

  P[, , 1] <- P0

  for(kk in 2:N){
    th_hat_min[, kk] <- A %*% th_hat[, kk-1]

    for(ii in 1:M){
      h_min[ii] <- sqrt(
        (th_hat_min[1,kk] - x[ii])^2
        + (th_hat_min[2,kk] - y[ii])^2)
      H[ii,] <- c((th_hat_min[1, kk] - x[ii])/h_min[ii],
                  (th_hat_min[2, kk] - y[ii])/h_min[ii], 0, 0)
    }

    # covariance estimation
    if(parameter$var.est == 1){
      sigma <- 1.483 * mean(abs((r_ges[,kk]-h_min)
                                - mean(r_ges[,kk]-h_min)))
      sigma2[kk] <- sigma^2
      R <- diag(sigma2[kk], M, M)
    }

    P_min[,,kk] <- A %*% P[,,kk-1] %*% t(A) + G %*% Q %*% t(G)
    K <- P_min[,,kk] %*% t(H) %*% inv(H %*% P_min[,,kk] %*% t(H) + R)

    th_hat[, kk] <- th_hat_min[,kk] + K %*% (r_ges[,kk] - h_min)
    P[,,kk] <- (diag(1, 4, 4) - K %*% H) %*% P_min[,,kk]
  }

  parameter$Rest <- sigma2
  parameter$K <- K

  return(list('th_hat' = th_hat,
              'P_min' = P_min,
              'P' = P,
              'parameter' = parameter))
}
