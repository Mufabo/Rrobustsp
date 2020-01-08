set.seed(1)

x <- rnorm(5)
P <- 1
#----
N <- length(x)
kap2 = 0.8724286 # kap=var(eta(randn(10000,1)))

# coarse grid search
phi_grid <- seq(-.99, .99, by = 0.05)

# finer grid via polynomial interpolation
fine_grid <- seq(-.99, .99, by = 0.001)

a_bip_sc <- numeric(length(phi_grid))
a_sc <- numeric(length(phi_grid))

# The following was introduced so as not to predict based on highly contaminated data in the
# first few samples.
x_tran <- x[1:min(10, floor(N / 2))]
sig_x_tran <- madn(x)
x_tran[abs(x_tran) > 3 * sig_x_tran] <-
  sign(abs(x_tran) > 3 * sig_x_tran) * 3 * sig_x_tran
x[1:min(10, floor(N / 2))] <- x_tran[1:min(10, floor(N / 2))]


bip_ar1_tau(x, N, phi_grid, fine_grid, kap2)
