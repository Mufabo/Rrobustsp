% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dependentData_armaEstBipM.R
\name{arma_est_bip_m}
\alias{arma_est_bip_m}
\title{arma_est_bip_m}
\usage{
arma_est_bip_m(x, p, q, beta_hat_s, a_sc_final)
}
\arguments{
\item{x:}{data (observations/measurements/signal)}

\item{p:}{autoregressive order}

\item{q:}{moving-average order}

\item{beta_hat_s:}{BIP S-estimate}

\item{a_sc_final:}{M scale estimate of residuals of BIP S-estimate}
}
\value{
ar_coeffs: vector of BIP-AR(p) MM-estimates

ma_coeffs: vector of BIP-MA(q) MM-estimates
}
\description{
The function  arma_est_bip_m(x,p,q) comuptes the BIP M-estimation step for BIP MM estimates of the
 ARMA model parameters. It can also be used as a stand-alone
 M-estimator.
}
\note{
File location: dependentData_armaEstBipM.R
}
\examples{
library(signal)
library(zeallot)
library(pracma)

N <- 500
a <- rnorm(N)
p <- 1
q <- 0
x <- signal::filter(1, c(1, -0.8), a)

beta_s <- arma_est_bip_s(x, p, q, tolX = 1e-8)
beta <- c(beta_s$ar_coeffs, beta_s$ma_coeffs)

arma_est_bip_m(x, p, q, beta, beta_s$inno_scale)
}
\references{
"Robust Statistics for Signal Processing"
  Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
  Cambridge University Press, 2018.

 "Bounded Influence Propagation \eqn{\tau}-Estimation: A New Robust Method for ARMA Model Estimation."
  Muma, M. and Zoubir, A.M.
  IEEE Transactions on Signal Processing, 65(7), 1712-1727, 2017.
}
