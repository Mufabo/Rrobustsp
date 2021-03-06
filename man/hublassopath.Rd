% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regression_hublassopath.R
\name{hublassopath}
\alias{hublassopath}
\title{hublassopath}
\usage{
hublassopath(y, X, c = NULL, L = 120, eps = 10^-3, intcpt = T,
  reltol = 1e-05, printitn = 0)
}
\arguments{
\item{L}{: Positive integer, the number of lambda values EN/Lasso uses.
Default is L=120.}

\item{reltol}{: Convergence threshold for IRWLS. Terminate when successive
estimates differ in L2 norm by a rel. amount less than reltol. default: 1e-05}

\item{y:}{Numeric data vector of size N x 1 (output, respones)}

\item{X:}{Numeric data matrix of size N x p. Each row represents one
observation, and each column represents one predictor (feature)
columns are standardized to unit length.}

\item{c:}{Threshold constant of Huber's loss function (optional;
otherwise use default value)}

\item{intcpt:}{Logical (true/false) flag to indicate if intercept is in the
regression mode. Default is true.}

\item{eps:}{Positive scalar, the ratio of the smallest to the
largest Lambda value in the grid. Default is eps = 10^-3.}

\item{printitn:}{print iteration number (default = 0, no printing)}
}
\value{
B    : Fitted M-Lasso regression coefficients, a p-by-(L+1) matrix,
           where p is the number of predictors (columns) in X, and L is
           the  number of Lambda values.

B0 : estimates values of intercepts

stats  : structure with following fields:
\itemize{
              \item Lambda = lambda parameters in ascending order
            \item sigma = estimates of the scale (a (L+1) x 1 vector)
            \item gBIC = generalized Bayesian information criterion (gBIC) value
                  for each lambda parameter on the grid.}
}
\description{
hublassopath computes the M-Lasso regularization path (over grid
of penalty parameter values) using Huber's loss function.
}
