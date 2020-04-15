#' create_environment_book
#'
#' The function computes generates data for tracking of a mobile user
#' equipment based on the values received by parameter. It is based on a function written by Ulrich Hammes,
#' Signal Processing Group, TU Darmstadt, February 2009
#'
#' @export
create_environment_book <- function(parameter, start, sigma_v){
  # to avoid no visible binding for global variable note
  rnorm <- NULL

  x <- nrow(parameter$BS)
  y <- ncol(parameter$BS)

  parameter$sigma.v <- sigma_v

  d <- matrix(0, parameter$M, parameter$N)

  # random force motion model
  xx <- matrix(0, parameter$dim*2, parameter$N)
  xx[,1] <- start

  for(ii in 1:(parameter$N-1)){
    xx[, ii+1] <- parameter$A %*% xx[,ii]
    + parameter$G %*% rnorm(2) * parameter$sigma.v[1,1]
  }

  parameter$truetrajectory <- xx

  for(jj in 1:parameter$M){
    for(ii in 1:parameter$N){
      d[jj, ii] <- sqrt((x[jj]-xx[1,ii])^2
                        + (y[jj]-xx[2,ii])^2)
    }
  }

  parameter$truedistances <- d

  # create noise
  nn <- matrix(0, parameter$M, parameter$N)
  index <- matrix(0, parameter$M, parameter$N)

  for(ii in 1:parameter$M){
    tmp <- markov_chain_book(parameter$MarkovChain[,,ii],
                             parameter$sigma.los,
                             parameter$sigma.nlos,
                             parameter$mu.nlos,
                             parameter$N)
    nn[ii,] <- tmp$y
    index[ii,] <- tmp$h
  }

  # alignment of variables
  parameter$start <- start
  parameter$truetrajectory <- xx
  parameter$truedistances <- d
  parameter$noise <- nn
  parameter$measureddistances <- d + nn

  ii <- parameter$numbermc
  parameter$thx <- parameter$truetrajectory[1,]
  parameter$thy <- parameter$truetrajectory[2,]
  parameter$thvx <- parameter$truetrajectory[3,]
  parameter$thvy <- parameter$truetrajectory[4,]

  return(parameter)
}
