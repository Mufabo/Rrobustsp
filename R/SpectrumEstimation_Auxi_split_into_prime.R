#'@export
split_into_prime <- function(y){
  N <- length(y)

  data('prime_numbers')

  primes <- prime_numbers$prime.numbers

  kk <- which(N < primes)[1]

  if(kk == 1){
    N_prime <- 1
  }
  else{
    if(N < tail(primes, 1)) N_prime <- primes[kk-1]
    if(N == tail(primes, 1)) N_prime <- tail(primes, 1)
    if(N > tail(primes, 1)){
      N_prime <- tail(primes, 1)
      warning('make longer list of prime numbers')
    }
  }
  if(N == N_prime) y_matrix <- y
  else y_matrix <- matrix(c(y[1:N_prime], y[(N-N_prime+1):N]), nrow = N_prime)
  return(y_matrix)
}
