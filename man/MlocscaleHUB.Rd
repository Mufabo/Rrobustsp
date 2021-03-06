% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Location_Scale.R
\name{MlocscaleHUB}
\alias{MlocscaleHUB}
\title{Huber's M-estimates of location and scale}
\usage{
MlocscaleHUB(y, c = NULL)
}
\arguments{
\item{y}{: real valued data vector of size N}

\item{c}{: tuning constant c>=0. Default is Null}
}
\value{
mu_hat : Huber's M-estimate of location

sigma_hat : Huber's M-estimate of scale

iter : integer. Number of iterations
}
\description{
Mlocscale computes Huber's M-estimates of location and scale.
}
\note{
File location : Location_Scale.R
}
\examples{
MlocscaleHUB(rnorm(5))

}
