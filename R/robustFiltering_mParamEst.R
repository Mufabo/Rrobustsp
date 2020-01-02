#' m_param_est
#'
#' The function computes an M-estimator of regression using the assymetric tanh score function.
#'
#'
#' @param recieve: is the received signal
#' @param Theta: contains the initial estimate
#' @param C: is the regression matrix of the model y = C*x + n
#' @param param.maxiters: maximal number of iterations
#' @param param.break: break condition
#'
#'
#' @return th2: M-estimate of regression
#' @return Theta:  M-estimates of regression for all iterations
#' @return kk: iteration index
#' @return residual: residuals given Theta
#'
#' @references
#'
#'   "Robust Statistics for Signal Processing"
#'   Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#'   Cambridge University Press, 2018.
#'
#'  "Robust Tracking and Geolocation for Wireless Networks in NLOS Environments."
#'   Hammes, U., Wolsztynski, E., and Zoubir, A.M.
#'   IEEE Journal on Selected Topics in Signal Processing, 3(5), 889-901, 2009.
#'
#'@examples
#'
#'@export
m_param_est <- function(receive, C, Theta, param){
  Pseudoinv <- MASS::ginv(C)
  resids <- c(receive - C %*% Theta)
  Theta <- cbind(numeric(length(Theta))
                 , Theta
                 , matrix(0, nrow = length(Theta), ncol = param$max.iters - 1)
                 )



  noisescale <- 1.483 * mean(abs(resids - mean(resids)))

  for(kk in 2:param$max.iters){
    thresh <- sum(abs((Theta[,kk] - Theta[,kk-1]) / Theta[,kk]))
    if(is.nan(thresh)) thresh <- Inf # number / 0 <- inf in matlab, but nan in R
    if(thresh < param$break.cond) break

    asym <- asymmetric_tanh(resids / noisescale
                            , param$c1
                            , param$c2
                            , param$x1)



    muu <- 1.25 * max(abs(asym$phi_point))

    # update param estimate
    Theta[,kk+1] <- Theta[,kk] +
      1 / muu * Pseudoinv %*% asym$phi * noisescale

    resids <- receive - C %*% Theta[,kk+1]

    noisescale <- 1.483 * mean(abs(resids - mean(resids)))
  }

  th2 <- Theta[,kk]

  return(list('th2' = th2,
              'Theta' = Theta,
              'kk' = kk-1,
              'residuals' = resids))
}
