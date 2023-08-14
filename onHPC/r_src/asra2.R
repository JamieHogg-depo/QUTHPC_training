## Model 1 ## -----------------------------------------------------------------

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
nI <- function(){list(Beta = rnorm(1),
                      alpha = rnorm(1),
                      sigma = runif(1, 0.1, 2))}
monitors <- c("alpha", "Beta", "sigma", "mu", "res")

## Fit and save the model ## ---------------------------------------------------

mod <- jf$runNimble(code = code, 
                           nD = nD,
                           nI = nI(),
                           nC = nC, 
                           monitors = monitors, 
                           niter = niter,
                           nburnin = nburnin, 
                           thin = thin, 
                           nchains = nchains)
						   
# Save mod object					   
saveRDS(mod, file = paste0(base_folder, "/outputs/", condition, "/", Rfile, "_mod.rds"))

# Summarise object 
fit <- jf$sumNimble(mod)

# remove mod object
rm(mod)

## Model fit and convergence ## -----------------------------------------------

# create list to save convergence plots
mf_cv_ll <- list()

# convergence diagnostic messages 
mf_cv_ll$conv_message <- fit$messages

# MCMC settings
mf_cv_ll$MCMCsettings <- fit$MCMCsettings

# Parameter summary
mf_cv_ll$param_summary <- fit$summary 

# trace plot
mf_cv_ll$traceplot <- mcmc_trace(as.list(fit$fit$samples), pars = c("alpha", "Beta", "sigma"))
# Add object NOT .png

# residual plot
res_draws <- jf$getSubsetDrawsNimble(as.matrix(fit$fit$samples), "res\\[")
fitted_draws <- jf$getSubsetDrawsNimble(as.matrix(fit$fit$samples), "mu\\[")

mf_cv_ll$resplot <-
data.frame(res =  apply(res_draws, 2, median),
		   fitt = apply(fitted_draws, 2, median)) %>%
ggplot(aes(y = res, x = fitt))+
  theme_bw()+
  geom_hline(yintercept = 0)+
  geom_point()+
  labs(y = "Residuals", 
       x = "Fitted values")
	   
# observed vs fitted
res_draws <- jf$getSubsetDrawsNimble(as.matrix(fit$fit$samples), "res\\[")
fitted_draws <- jf$getSubsetDrawsNimble(as.matrix(fit$fit$samples), "mu\\[")

mf_cv_ll$fittedplot <- 
data.frame(y = data$y,
		   fitt = apply(fitted_draws, 2, median)) %>%	
ggplot(aes(y = y, x = fitt))+
  theme_bw()+
  geom_abline()+
  geom_point()+
  labs(y = "Observed", 
       x = "Fitted values")
	   
# save convergence plots
saveRDS(mf_cv_ll, file = paste0(base_folder, "/outputs/", condition, "/", Rfile, "_convergence.rds"))
						   
## Summarise draws ## ---------------------------------------------------------

results1 <- jf$getResultsData(fitted_draws, 
                              model = paste0(condition, "_", sex), 
                              metric = "standard")
saveRDS(results1, file = paste0(base_folder, "/outputs/", condition, "/", Rfile, "_results1.rds"))
							
results2 <- jf$getResultsData(exp(fitted_draws), 
                              model = paste0(condition, "_", sex), 
                              metric = "exponentiated")
saveRDS(results2, file = paste0(base_folder, "/outputs/", condition, "/", Rfile, "_results2.rds"))

## END SCRIPT ## --------------------------------------------------------------



