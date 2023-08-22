
# Progress bar for this: https://masonfidino.com/nimble_parallel_pb/

this_cluster <- makeCluster(4)

code <- nimbleCode({
  for (i in 1:N) {
    # linear predictor
    mu[i] <- alpha + Beta * x[i]
    # likelihood
    y[i] ~ dnorm(mu[i], sd = sigma)
  }
  
  # Priors
  alpha ~ dnorm(0, sd = 100)
  Beta ~ dnorm(0, sd = 100)
  sigma ~ dgamma(2, 0.5)
  
  # residuals
  res[1:N] <- y[1:N] - mu[1:N]
})

run_MCMC_allcode <- function(seed, data, code){
  
  library(nimble)
  source('C:/r_proj/QUTHPC_training/onHPC/r_src/functions.R')
  
  ## Fit and save the model ## ---------------------------------------------------
  
  mod <- jf$runNimble(code = code, 
                      nD = list(y = data$y,
                                x = data$x),
                      nI = function(){list(Beta = rnorm(1),
                                           alpha = rnorm(1),
                                           sigma = runif(1, 0.1, 2))},
                      nC = list(N = nrow(data)), 
                      monitors = c("alpha", "Beta", "sigma", "mu", "res"), 
                      niter = 100,
                      nburnin = 50, 
                      thin = 1, 
                      nchains = 1)
  
  return(mod)
}

# Run parallel
chain_output <- parLapply(cl = this_cluster, X = 1:4,
                          fun = run_MCMC_allcode,
                          data = data, code = code)
stopCluster(this_cluster)




