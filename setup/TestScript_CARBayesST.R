# Set library
	.libPaths("r_lib")

# Install package
	install.packages("CARBayesST")

# Load CARBayesST
  library(CARBayesST)

# simulate data
  N <- 10
  T <- 10
  pop <- runif(N*T, 300, 600)
  rate <- runif(N*T, 0.01, 0.1)
  y <- rpois(N*T, rate*pop)
  W <- matrix(0, ncol = N, nrow = N)
  W[1,2] <- 1
  W[2,c(1,3)] <- 1
  W[3,c(2,4)] <- 1
  W[4,c(3,5)] <- 1
  W[5,c(4,6)] <- 1
  W[6,c(5,7)] <- 1
  W[7,c(6,8)] <- 1
  W[8,c(7,9)] <- 1
  W[9,c(8,10)] <- 1
  W[10,9] <- 1

# Fit the Bayesian model
  fit <- ST.CARanova(y ~ offset(log(pop))+1, 
                     family = "poisson",
                     W = W, 
                     burnin = 2000,
                     n.sample = 4000)

# Summarize results
  fit$summary.results
