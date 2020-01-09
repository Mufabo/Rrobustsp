#'order_wk
#'
#'@examples
#'
#'library(signal)
#'library(Rrobustsp)
#'library(zeallot)
#'
#'y <- rnorm(5)
#'
#'c(wr, PSD) <- order_wk(y)
#'
#'@note
#'
#'Location: .../Rrobuststp/R/SpectrumEstimation_Auxi_order_wk.r
#'
#'@export
order_wk <- function(y) {
  yt <- y - median(y)


  N <- length(y)
  K <- (N - 1) / 2

  ARM <- numeric(K)
  BRM <- numeric(K)
  PSD <- numeric(K)

  w <- numeric(K)

  for (k in 1:K) {
    w[k] <- 2 * pi * k / N

    Apuv <- matrix(0, N, N)
    Bpuv <- matrix(0, N, N)

    Ap <- numeric(K)
    Bp <- numeric(K)
      for (u in 0:(N - 1)) {
        for (v in 0:(N - 1)) {
          if (u != v) {
            denom <- sin(w[k] * (v - u))
            Apuv[u+1, v+1] <- ( yt[u+1] * sin(w[k] * v) - yt[v+1] * sin(w[k] * u)) / denom
            Bpuv[u+1, v+1] <- ( yt[v+1] * cos(w[k] * u) - yt[u+1] * cos(w[k] * v)) / denom
          }
        }
        Ap[k] <- median(apply(Apuv, 1, median))
        Bp[k] <- median(apply(Bpuv, 1, median))
        PSD[k] <- Ap[k]^2 + Bp[k]^2
      }
  }

  N_win <- 3
  Smoothing_Win <- signal::hamming(N_win)
  PSD_smooth <- signal::conv(PSD, Smoothing_Win)
  PSD_smooth <- head(tail(PSD_smooth, -N_win), -N_win)

  c(PSD_sort, order) %<-% sort(PSD, decreasing = T, index.return = T)
  wr <- 2 * pi * order / N

  return(list(wr, PSD))
}
