
repeated_median_filter <- function(x){
  N <- length(x)
  xFRM <- numeric(N)

  x_split <- split_into_prime(x)

  for(ii in 1:ncol(x_split)){
    x_part <- x_split[,ii]

    wr <- order_wk(x_part)

    N_prime <- length(x_part)

    t <- 0:(N_prime - 1)
    K <- (N_prime - 1) / 2

    ARM <- numeric(K)
    BRM <- numeric(K)

    M <- 2

    xm <- median(x_part)
    xt <- x_part - xm

    for(m in 1:M){
      for(k in 1:K){
        for(u in 0:(N_prime - 1)){
          for(v in 0:(N_prime - 1)){
            if(u != v){
              denom <- sin(w[k] * (v - u))
              Apuv[u+1, v+1] <- ( yt[u+1] * sin(w[k] * v) - yt[v+1] * sin(w[k] * u)) / denom
              Bpuv[u+1, v+1] <- ( yt[v+1] * cos(w[k] * u) - yt[u+1] * cos(w[k] * v)) / denom
            }
          }
          A[k] <- median(apply(Apuv, 1, median))
          B[k] <- median(apply(Bpuv, 1, median))
          xt <- xt - A[k] * cos(wr[k] * t) - B[k] * sin(wr[k] * t)
          ARM[k] <- ARM[k] + A[k]
          BRM[k] <- BRM[k] + B[k]
        }
      }
    }

    if(ii == 1){
      xFRM1 <- xm
      for(k in 1:K){
        sumAB <- ARM[k] * cos(wr[k]*t) + BRM[k]*sin(wr[k]*t)
        xFRM1 <- xFRM1 + sumAB
      }
    }

    if(ii == 2){
      xFRM2 <- xm
      for(k in 1:K){
        sumAB <- ARM[k] * cos(wr[k]*t) + BRM[k]*sin(wr[k]*t)
        xFRM2 <- xFRM2 + sumAB
      }
    }
  }

  if(ii==1) xFRM <- xFRM1
  else{
    xFRM[1:(N - N_prime)] <- xFRM1[1:(N - N_prime)]
    xFRM[(N-N_Prime+1):N_prime] <- 0.5 * (xFRM1[(N-N_Prime+1):N_prime] + xFRM2[(N-N_Prime+1):N_prime])
    xFRM[(N_prime+1):N] <- xFRM2[(N_prime+1):N]
  }

  return(list('xFRM' = xFRM, 'ARM' = ARM, 'BRM' = BRM))
}
