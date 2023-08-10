## Explore unit-level survey data ## -------------------------------------------

# load data
load(file = paste0(base_loc, "HWSS_smoke_unit_data.Rdata"))
df <- HWSS_smoke_unit_data$df; poststrata <- HWSS_smoke_unit_data$poststrata
rm(HWSS_smoke_unit_data, hr, hd, sa2)

# Set parameters for MCMC
MCMCsettings <- list(print_samplers = F,
						optimBeta = T,
						beta_name = paste0("B_qr[", 1:16, "]"),
						sampler_name = "ess",
						# increase the frequency of adaption during burnin
						adaptInterval = 20)

## Derive direct estimates for years and areas ## ------------------------------

# get population counts by area and year
  pop <- poststrata %>% 
    group_by(MT_id) %>% 
    summarise(N = sum(N), .groups = "drop")
  
# get the direct estimates
  d <- jf$jDirect(df$y, 
                  df$MT_id,
                  sweight = df$w_ss,
                  domsize = pop,
                  area_name = "MT_id")
  
## Use Non-Bayesian MrP --------------------------------------------------------

# fit a weighted model
  fit_w <- glm(y ~ 1,
               family = binomial,
               data = df,
               weights = w_ss) # include weights
  s_fit_w <- summary(fit_w)
  
  # use the fit to get good starting values for the fixed effects in MCMC
  alpha_prior <- list(mean = s_fit_w$coefficients[1,1], se = s_fit_w$coefficients[1,2])
  
## Setup for Nimble modelling ## -----------------------------------------------
  
  # Load new nimble functions
    # vectorized pseudo-likelihood functions
  dwbern_v <- nimbleFunction(
    run = function(x = double(1), 
                   w = double(1), 
                   p = double(1), 
                   log = integer(0, default = 0)) {
      returnType(double(0))
      #logProb <- sum(w * dbinom(x, size = rep(1, length(x)), prob = p, log = TRUE))
      logProb <- sum(w * (x * log(p) + (1-x) * log(1-p)))
      if(log) return(logProb)
      else return(exp(logProb)) 
    })
  # unfortunately we also have to load a function for the the random function
    # NOTE: This function should NEVER be used! 
  rwbern_v <- nimbleFunction(
    run = function(n = integer(0), 
                   w = double(1), 
                   p = double(1)) {
      returnType(double(1))
      return(w * rbinom(n, size = rep(1, length(p)), prob = p))
    })
  
# derive the spatial weight matrix
  map <- lga$map %>% 
    # remove empty geometries
    filter(!st_is_empty(.)) %>% 
    # include the concordance id, M_id
    right_join(.,(poststrata %>% group_by(M_id, LGA_NAME16) %>% tally() %>% dplyr::select(-n)), 
               by = "LGA_NAME16") %>% 
    # ensure the order is correct
    arrange(M_id) %>% 
    # ensure the result is s simple features (sf) object
    st_as_sf()
  # get the nb list
  nb <- poly2nb(map)
  # get the weight matrix
  M_W <- nb2mat(nb, style="B")
  # derive the necessary scaling factor for W
  scaling <- jf$getBYM2scale(M_W)
  # get the needed elements from W
  M_bugs <- as.carAdjacency(M_W)
  # get the W elements for the temporal W matrix
  T_bugs <- jf$tm.adjacency(T = 9, 1)

# Construct design matrix and take QR decomposition
  # SURVEY DATA
  X_s <- model.matrix(~agegroup*sex + RA_Name + IRSD_5 + age1840 + age4065 + female_prop, data = df)
  QR_s <- jf$getQRDecomp(X_s)
  # POSTSTRATA DATA
  QR_ps <- jf$getQRDecomp(model.matrix(~agegroup*sex + RA_Name + IRSD_5 + age1840 + age4065 + female_prop, data = poststrata),
                          center = colMeans(X_s))

## Bayesian fit ## -------------------------------------------------------------

# Nimble/BUGS model
code <- nimbleCode({
  # vectorized pseudo-likelihood for bernoulli
  y[1:N] ~ dwbern_v(w = w[1:N], p = pr[1:N])
  
  for(i in 1:N){
    # linear predictor
    logit(pr[i]) <- alpha + inprod(B_qr[1:q], Q_ast[i,]) + s[M_id[i]] + t[T_id[i]]
  }
  
  # Spatial: BYM2
  for(i in 1:M){
    s[i] <- sigma_s * (Corr[i] + UCorr[i])
    # structured component
    Corr[i] <- Z_s[i] * sqrt(rho/scale)
    # unstructured component
    UCorr[i] <- Z_v[i] * sqrt((1 - rho))
    # standard normal
    Z_v[i] ~ dnorm(0, 1)
  } 
  
  # Spatial: ICAR prior
  Z_s[1:M] ~ dcar_normal(adj[1:L_s], weights[1:L_s], num[1:M], 1, zero_mean = 1)

  # Temporal: RW1 prior
  Z_t[1:T] ~ dcar_normal(T_adj[1:L_t], T_weights[1:L_t], T_num[1:T], 1, zero_mean = 1)
  t[1:T] <- sigma_t * Z_t[1:T]
  
  # Other priors
  for(i in 1:q){
    B_qr[i] ~ dnorm(0,2)
  }
  alpha ~ dflat()
  sigma_s ~ dgamma(2, 0.5)
  sigma_t ~ dgamma(2, 0.5)
  rho ~ dunif(0,1)
  
  # recreate true coefficient values
  B[1:q] <- R_ast_inverse[1:q,1:q] %*% B_qr[1:q]
})

# data and constant lists
nD <- list(y = df$y,
           w = df$w_ss,
           Q_ast = QR_s$QR$Q_ast,
           R_ast_inverse = QR_s$QR$R_ast_inverse)
nC <- list(q = QR_s$p,
           N = nrow(df),
           M = length(unique(df$M_id)),
           M_id = df$M_id,
           T = length(unique(df$T_id)),
           T_id = df$T_id,
           L_s = length(M_bugs$adj),
           L_t = length(T_bugs$T_adj),
           scale = scaling)
# add adjacency objects to constant list
nC <- c(nC, M_bugs, T_bugs)
# initial values
nI <- function(){list(Z_s = jf$vectorSumtoZero(nC$M),
                      Z_t = jf$vectorSumtoZero(nC$T),
                      Z_v = rnorm(nC$M),
                      B_qr = rnorm(nC$q),
                      sigma_s = runif(1, 0.2, 2),
                      sigma_t = runif(1, 0.2, 2),
                      rho = runif(1, 0,1),
                      # use specific initial values from MLE
                      alpha = rnorm(1, alpha_prior$mean, alpha_prior$se))}
monitors <- c("sigma_s", "sigma_t", "B", "B_qr", "t", "alpha", "rho", "s", "Corr")

## Fit the weighted MRP model ## ----------------------------------------------
mrp_fit <- jf$SampleNimble(code = code, 
                           nD = nD,
                           nI = nI(),
                           nC = nC, 
                           monitors = monitors,
                           niter = niter,
                           nburnin = nburnin, 
                           thin = thin, 
                           nchains = nchains, 
                           print_samplers = MCMCsettings$print_samplers,
                           optimBeta = MCMCsettings$optimBeta,
                           beta_name = MCMCsettings$beta_name,
                           sampler_name = MCMCsettings$sampler_name,
                           # increase the frequency of adaption during burnin
                           adaptInterval = MCMCsettings$adaptInterval)
						   
## Get posterior draws - MRP ## -----------------------------------------------------

# Extract posterior draws for ALL parameters
fit_draws <- as.matrix(mrp_fit$fit$samples)
# get draws for specific parameters
  alpha <- fit_draws[,which(str_detect(attr(fit_draws, "dimnames")[[2]], "alpha"))]
  beta <- fit_draws[,which(str_detect(attr(fit_draws, "dimnames")[[2]], "B_qr\\["))]
  # both spatial and temporal random effects - iterations by MxT
  lambda <- cbind(fit_draws[,which(str_detect(attr(fit_draws, "dimnames")[[2]], "s\\["))], 
                  fit_draws[,which(str_detect(attr(fit_draws, "dimnames")[[2]], "t\\["))])

# construct the sparse Z matrix for poststrata
Z <- Matrix(cbind(
  model.matrix(~as.factor(M_id) - 1, data = poststrata),
  model.matrix(~as.factor(T_id) - 1, data = poststrata)
), sparse = T)

# prepare matrix of draws and run for loop
  # empty matrix of iterations by length(poststrata)
  counts <- matrix(NA, nrow = length(alpha), ncol = nrow(poststrata))
  # add a progress bar
  pb <- txtProgressBar(min = 0, max = length(alpha), style = 3)
  for(i in 1:length(alpha)){
    # right side returns a vector of length nrow(poststrata)
    counts[i,] <- poststrata$N * as.numeric(jf$jinvlogit(alpha[i] + QR_ps$QR$Q_ast %*% beta[i,] + Z %*% lambda[i,]))
    setTxtProgressBar(pb, i)
  }
  close(pb) # remove the progress bar
  
# get posterior draws by area and year
  temp <- poststrata %>% group_by(MT_id) %>% summarise(N = sum(N))
  mrp_prev_draws <- t(apply(counts, 1, FUN = function(x){aggregate(x, list(poststrata$MT_id), sum)[,2]/temp$N}))
    # ordered according to MT_id
    # iterations x MT
  rm(temp) # remove temp
  
# get results
  # user made function that returns point estimates and intervals
  # we also add difference in posterior probabilities too with addDPP
  mrp_results <- jf$getResultsData(mrp_prev_draws, 
                    model = "MRP", 
                    metric = "Prevalence",
                    prefix = "mrp_",
                    other_data = poststrata %>% group_by(MT_id,year, LGA_NAME16) %>% 
                      tally() %>% dplyr::select(-n),
                    addDPP = T, 
                    null_value = 0.09450604)

## SAVE environment ## --------------------------------------------------------

out <- list(fit = mrp_fit, 
            prev_draws = mrp_prev_draws,
            results = mrp_results)
saveRDS(out, file = paste0("WADOH/outputs/", cur_date, "/r/", Rfile, ".rds"))



