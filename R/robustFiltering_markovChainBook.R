markov_chain_book <- function(MM, sigma_los, sigma_nlos, mu_nlos, N){

  # transition probabilities
  pt_losnlos <- MM[1,2]
  pt_nloslos <- MM[2,1]


  h <- numeric(length = N+100)
  y <- numeric(length = N+100)

  pp <- runif(1)

  # initialisation
  if(pp < 0.5){
    state <- 0 # LOS
  } else {
    state <- 1 # NLOS
  }

  for(m in 1:(N + 100)){
    if(state == 0){
      # LOS
      p <- runif(1)
      if(m > 100){
        # ensures that the Markov Chain is in a steady state (p_(k) = p_init*A^k), if k large, p(k) is steady, where A is the Markov transition matrix
        y[m] <- sigma_los * rnorm(1)
        h[m] <- 0
      }
      if(p < pt_losnlos) state <- 1
    } else {
      # NLOS
      q <- runif(1)
      if(m > 100){
        y[m] <- sigma_los * rnorm(1) + mu_nlos + sigma_los * rnorm(1)
        h[m] <- 1
      }
      if(q < pt_nloslos) state <- 0
    } # end else
  } # end for

  # the first 100 samples are discarded to ensure that the Markov Chain is in
  # the steady state

  y <- y[101:length(y)]
  h <- h[101:length(h)]

  return(list('y' = y, 'h' = h))
} # end function
