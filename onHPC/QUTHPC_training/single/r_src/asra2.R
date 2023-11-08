## Model 2 - parallel ## ------------------------------------------------------

library(doParallel)

# load data
data <- readRDS(file = paste0(base_folder, "/data/", model_spec, "_", condition, ".rds"))

# filter by sex
data <- data %>% 
  filter(sex == sex)

## Setup model ## -------------------------------------------------------------

code <- nimbleCode({
	for (i in 1:N) {
	  # linear predictor
	  mu[i] <- alpha + Beta * x[i]
	  # likelihood
	  y[i] ~ dnorm(mu[i], sd = sigma)
	}

	# Priors
	alpha ~ dnorm(0, sd = 1)
	Beta ~ dnorm(0, sd = 1)
	sigma ~ dgamma(2, 0.5)
	
	# residuals
	res[1:N] <- y[1:N] - mu[1:N]
})

nD <- list(y = data$y,
           x = data$x)
nC <- list(N = nrow(data))
monitors <- c("alpha", "Beta", "sigma", "mu", "res")

## Fit and save the model ## ---------------------------------------------------

# open cluster
this_cluster <- makeCluster(4)

# define parallel function
run_MCMC_allcode <- function(seed, base_folder, code, nD, nC, monitors, niter, nburnin, thin){
  
  .libPaths('r_lib')
  library(nimble)
  source(paste0(base_folder, "/r_src/functions.R"))
  
  mod <- jf$runNimble(code = code, 
                      nD = nD,
                      nI = function(){list(Beta = rnorm(1),
                      alpha = rnorm(1),
                      sigma = runif(1, 0.1, 2))},
                      nC = nC, 
                      monitors = monitors, 
                      niter = niter,
                      nburnin = nburnin, 
                      thin = thin, 
                      nchains = 1)
  
  return(mod)
}

# Run model in parallel
mod <- parLapply(cl = this_cluster, X = 1:4,
fun = run_MCMC_allcode,
base_folder=base_folder,
code = code, 
nD = nD,
nC = nC, 
monitors = monitors, 
niter = niter,
nburnin = nburnin, 
thin = thin)

# close cluster
stopCluster(this_cluster)
						   
# Save mod object					   
saveRDS(mod, file = paste0(loop_output_file, "_mod.rds"))

## END SCRIPT ## --------------------------------------------------------------



