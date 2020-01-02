parameter <- list()

parameter$N <- 3000 # length of simulated trajector
parameter$BS <- matrix(c(2000, 7000, 12000, 7000, 7000, 12000, 7000, 2000, 7000, 7000), 5, 2, byrow = T)
parameter$M <- nrow(parameter$BS) # number of Base stations
parameter$start <- c(4300, 4300, 2, 2) # mean starting position
parameter$pnlos <- c(0,0,0,0,0.25)
parameter$initial_sigma <- c(50, 50, 6, 6) # standard deviation for the for the initial state
parameter$mc <- 50 # number of monte carlo runs
parameter$discardN <- 100
parameter$grid <- 1
parameter$Ts <- 0.2
parameter$A <- matrix(data = c(1.,0,0,0
                               , 0,1,0,0
                               , parameter$Ts,0,1,0,
                               0, parameter$Ts,0,1), 4, 4)
parameter$G <- rbind(parameter$Ts^2 / 2 * diag(1, 2, 2)
                     , parameter$Ts * diag(1,2,2))


parameter$dim <- ncol(parameter$BS)
parameter$numberbs <- 'variable'
parameter$noisemodel <- 'GMM1'
parameter$motionmodel <- 'random-force'
parameter$noisestructure <- 'Markov'
parameter$figure <- 1
parameter$plot <- 'mse'
parameter$MarkovChain <- to.tensor(0,c(U=2
                                       ,V=2
                                       ,W=parameter$M))







#save.image(file = '~/Rrobustsp/data/set_parameters_book.RData')
