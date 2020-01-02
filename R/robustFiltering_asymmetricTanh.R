#' asymmetric_tanh
#'
#' Computes the smoothed asymmetric tanh score
#' and its derivative.
#'
#' @param c1 : first clipping point
#' @param c2 : second clipping point
#' @param x1 : smoothing parameter to make the score function continuous
#' @param sig : the signal
#'
#' @return phi : asymmetric tanh score function
#' @return phi_point : derivative of asymmetric score function
#'
#' @examples
#'
#' asymmetric_tanh(c(1, 3, 4, 5, 62, 346, 23562, 65, 264, 21362, 90), 6, 87, 4)
#' @references
#' "Robust Statistics for Signal Processing"
#'  Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#'  Cambridge University Press, 2018.
#'
#'  "Robust Tracking and Geolocation for Wireless Networks in NLOS Environments."
#'   Hammes, U., Wolsztynski, E., and Zoubir, A.M.
#'   IEEE Journal on Selected Topics in Signal Processing, 3(5), 889-901, 2009.
#' @export
asymmetric_tanh <- function(sig, c1, c2, x1){
  phi <- numeric(length(sig))

  phi[abs(sig) <= c1] <- sig[abs(sig) < c1]
  phi[abs(sig) >  c1] <- x1 * mat_sign(sig[abs(sig) > c1]) *
    tanh(x1 * 0.5 * (c2 - abs(sig[abs(sig) > c1])))
  phi[abs(sig) > c2] <- numeric(length = length(sig[abs(sig) > c2]))

  phi_point <- numeric(length(sig))
  phi_point[abs(sig) <= c1] <- 1
  phi_point[abs(sig) > c1] <- -0.5 * x1^2 /
    (cosh(0.5 * x1 * (c2 - abs(sig[abs(sig) > c1]))))^2
  phi_point[abs(sig) > c2] <- numeric(length = length(sig[abs(sig) > c2]))

  return(list('phi' = phi, 'phi_point' = phi_point))

}
