# Set library
.libPaths("r_lib")

# Install package
install.packages("nimble", repos = "http://r-nimble.org", type = "source")

# Load Nimble
library(nimble)

# simulate data
N <- 100
y <- rnorm(N, 1, 0.2)

# write the code
code <- nimbleCode({
  # Model
  for (i in 1:N) {
    y[i] ~ dnorm(mu, sd = sigma)
  }
  
  # Priors
  mu ~ dnorm(0,1)
  sigma ~ dgamma(2,0.5)
})

# define the constants
cons <- list(N = N)

# define the data
data <- list(y = y)

# initial values
inFun <- function(){
  list(sigma = runif(1, 0.1, 0.4), mu = rnorm(1))
}

# Define the model object
fit <- nimbleMCMC(code = code,
                  data = data,
                  inits = inFun(),
                  constants = cons,
                  nchains = 4, 
                  nburnin = 2000,
                  summary = T,
                  WAIC = T)

# Summarize
fit$summary$all.chains

