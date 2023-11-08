## Functions ## ----------------------------------------------------------------

## Create master list
jf <- list()

## ----------------------------------------------------------------------------
jf$givemenicetitles <- function(input){
	wrapped <- paste0("## ", input, " ## ")
	nu_lines <- 80 - str_length(wrapped)
	message(paste0(wrapped, paste0(rep("-", nu_lines), collapse = "")))
}

## ----------------------------------------------------------------------------
jf$lu <- function(x){
	length(unique(x))
}

## ----------------------------------------------------------------------------
jf$jggsave <- function(filename, path, dims = c(3000, 2300)){
    ggsave(filename = filename, path = path,
           dpi = 300,
           width = dims[1],
           height = dims[2],
           units = "px")
}

## ----------------------------------------------------------------------------
#' @description TT stands for tidyTable
jf$TT <- function(.data, ..., .include_n = T){
	if(.include_n){
		.data %>% 
			group_by(...) %>% 
			tally() %>% ungroup()
	}else{
		.data %>% 
			group_by(...) %>% 
			tally() %>% ungroup() %>%
			dplyr::select(-n)
	}
}

## ---------------------------------------------------------------------------
#' @description get size of all objects in R session
jf$ObjectSizeCheck <- function(){
	data.frame('object' = ls()) %>% 
	  dplyr::mutate(size_unit = object %>%sapply(. %>% get() %>% object.size %>% format(., unit = 'auto')),
					size = as.numeric(sapply(strsplit(size_unit, split = ' '), FUN = function(x) x[1])),
					unit = factor(sapply(strsplit(size_unit, split = ' '), FUN = function(x) x[2]), levels = c('Gb', 'Mb', 'Kb', 'bytes'))) %>% 
	  dplyr::arrange(unit, dplyr::desc(size)) %>% 
	  dplyr::select(-size_unit) %>% 
	  relocate(unit, size)
}

## ----------------------------------------------------------------------------
#' @description requires the xtable packages
#' @param digits
#' @returns latex code for a table. Function can be used in pipe. 
# Any number of variables can be supplied to the `...`
jf$getLatexTable <- function(.data, digits = 2, ...){
  .data %>% 
    dplyr::select(...) %>% 
    relocate(...) %>% 
    head(n = 11) %>% 
    xtable(., digits = digits)
}

## ----------------------------------------------------------------------------
jf$Mat2ListbyRows <- function(mat){
	if(!is.matrix(mat)){
		split(mat, f = seq(length(mat)))
	}else{
		split(mat, f = seq(nrow(mat)))
	}
}

## ----------------------------------------------------------------------------
#' @param p vector of probabilities to convert to unconstrained space
#' @upper numeric scalar (defaults to 1). Useful when using a adjusted logit transformation
jf$jlogit <- function(p, upper = 1){
  log(p/(upper-p))
}

## ----------------------------------------------------------------------------
#' @param eta vector of unconstrained variables to convert between 0 and upper
#' @upper numeric scalar (defaults to 1). Useful when using a adjusted logit transformation
jf$jinvlogit <- function(eta, upper = 1){
  (upper * exp(eta))/(1+exp(eta))
}

## ----------------------------------------------------------------------------
jf$jUnique <- function(x){
  r <- unique(x)
  ifelse(is.na(r), "NA", r)
}

## ----------------------------------------------------------------------------
#' returns data.frame where each row is a variable with the number of unique elements
# and what the unique elements are. The function does NOT show variables that are
# constant (i.e. `length(unique(x)) == 1`)
jf$getUnique <- function(.data){
  .data %>% 
    map(~str_c(jf$jUnique(.x), collapse = ",")) %>% 
    bind_rows() %>% 
    pivot_longer(cols = everything()) %>% 
    mutate(n = str_count(value, ",") + 1) %>% 
    relocate(name, n) %>% 
    filter(n > 1)
}

## ----------------------------------------------------------------------------
#' @param name non-character for the variable that will be added
#' @param concor logical (defaults to F). If using the function to create a separate
# concordance dataset set to T, otherwise the function adds the column to the existing data.
#' @param all_combs logical (defaults to T). Only applicable when concor = T. Setting to F
# only generates an ID for the combinations in the given data.  
# the dots are for the group variable
jf$addGroupID <- function(.data, name, ..., concor = F, all_combs = T){
  if(concor){
	if(all_combs){
		.data %>% 
			dplyr::select(...) %>% 
			filter(!duplicated(.)) %>% 
			complete(...) %>%
			group_by(...) %>% 
			summarise("{{name}}" := cur_group_id(),
						.groups = "drop") %>%
			filter(complete.cases(.))
	}else{
	.data %>% 
      group_by(...) %>% 
      summarise("{{name}}" := cur_group_id(),
                .groups = "drop")
	}
  }else{
    .data %>% 
      group_by(...) %>% 
      mutate("{{name}}" := cur_group_id()) %>% 
      ungroup() 
  }
}

## ----------------------------------------------------------------------------
#' @description filters by the first unique combination of the grouping factors given
# Requires the function addGroupID
jf$checkGroupCombo <- function(.data, ..., index = 1){
	.data %>% 
	  jf$addGroupID(.temp, ...) %>% 
		  filter(.temp == index) %>% 
		  dplyr::select(-.temp)
}

## ----------------------------------------------------------------------------
jf$all.na <- function(x){
  length(which(is.na(x))) == length(x)
}

## ----------------------------------------------------------------------------
#' @param U number of random draws per observation
#' @param mu vector of length n_obs
#' @param var vector of length n_obs
#' @param pop vector of length n_obs
#' @return list of two: 
# 		y (n_obs x U) -> non-rounded counts
#		y_tilde matrix (n_obs x U) -> ceiling() counts
# 		offset_N_tilde matrix (n_obs x U)
jf$rbetaMP <- function(U, mu, var, pop) {

	n_obs <- length(mu)

	# get parameters of beta
	phi = ((mu * (1-mu))/var) - 1
	stopifnot("Check constraints of mu and var..." = all(phi >= 0))
	
	alpha = mu * phi
	beta <- (1-mu) * phi
	# alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
	# beta <- alpha * (1 / mu - 1)

	# get random draws
	draws <- rbeta(U*n_obs, shape1 = alpha, shape2 = beta)

	# get counts
	y <- draws * pop
	y_tilde <- ceiling(y)
	c <- ifelse(y_tilde/y == 0, 1, y_tilde/y) 

	# return the matrix of y_tilde
	return(list(y = matrix(y, byrow = F, ncol = U, nrow = n_obs),
				y_tilde = matrix(y_tilde, byrow = F, ncol = U, nrow = n_obs),
				N_tilde = matrix(c*pop, byrow = F, ncol = U, nrow = n_obs))
			)
}

## ----------------------------------------------------------------------------
#' @param n number of draws to return
#' @returns a vector of length n where the sum is equal to zero
jf$vectorSumtoZero <- function(n){
  init <- numeric(n)
  init[1:n-1] <- rnorm(n-1)
  init[n] <- 0 - sum(init[1:n-1])
  return(init)
}

# ## ----------------------------------------------------------------------------
# # Load age brackets
# age_brc <- data.frame(start = seq(0,80,5),
					  # end = seq(4,84, 5)) %>% 
					  # mutate(age = paste0(start, "_", end))
# age_labs <- data.frame(name = 1:19,
					 # age = c(age_brc$age, "85+", "Total"))
# rm(age_brc)

# ## -----------------------------------------------------------------------------
# hr_labs <- data.frame(hr_id = 1:10,
                      # hr_name = c("East MHS",
                                  # "North MHS",
                                  # "South MHS",
                                  # "Goldfields HR",
                                  # "Great Southern HR",
                                  # "Kimberley HR",
                                  # "Midwest HR",
                                  # "Pilbara HR",
                                  # "South West HR",
                                  # "Wheatbelt HR"))

## -----------------------------------------------------------------------------
#' @param x vector
#' @param range sequence of integers. Best to use 0:4 syntax
#' @param replaceWhat numeric or character describing when values 
# should be randomly replaced
jf$randomInts <- function(x, range, replaceWhat){
  
  temp <- x
  
  for(i in 1:length(x)){
    temp[i] <- ifelse(x[i] == replaceWhat, sample(range, 1), x[i])
  }
  return(temp)
}

## -----------------------------------------------------------------------------
#' @describe Will remove any columns where values are constant. 
# Best to use in a pipe setting `data %>% dropConstantCols()`
#' @param verbose logical (defaults to F). Setting to true
# prints the dropped variables to the console. 
jf$dropConstantCols <- function(.data, verbose = F){
  idss <- .data %>% 
    map(~str_c(unique(.x), collapse = ",")) %>% 
    bind_rows() %>% 
    pivot_longer(cols = everything()) %>% 
    mutate(n = str_count(value, ",") + 1) %>% 
    relocate(name, n) %>% 
    filter(n == 1)
  to_remove <- which(colnames(.data) %in% idss$name)
  
  # print message
  if(verbose){
	message("Dropped ", paste(idss$name, collapse = ", "))
  }
  
  # return data
  .data %>% 
    dplyr::select(-all_of(to_remove))
}

## -----------------------------------------------------------------------------
#' @param x vector of SIRs used as fill in geom_sf
#' @param extreme numeric specifying the largest SIR for which data is collapsed
#' @param nu_cats integer for the total number of labelled categories (must be odd)
jf$addJessColorScaleSIR <- function(x, extreme = 2, nu_cats = 7){
  
  SIR.clip <- x
  nu_cats_below <- (nu_cats - 1)/2
  
  # Trim extreme values
  cut.offs <- c(1/extreme, extreme)
  #set extreme values to bounds
  SIR.clip[which(SIR.clip < cut.offs[1], arr.ind = TRUE)] <- cut.offs[1]
  SIR.clip[which(SIR.clip > cut.offs[2], arr.ind = TRUE)] <- cut.offs[2]
  
  # Fill colours
  Fill.colours <- c("#2C7BB6", "#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C", "#D7191C")
  End <- log2(extreme + 0.1)
  
  # Create breaks
  Breaks.fill <- c(seq(1/extreme, 1, length = nu_cats_below)[-nu_cats_below],
                   1, 
                   seq(1, extreme, length = nu_cats_below)[-1])
  Fill.values <- c(-End, log2(Breaks.fill), End)
  Fill.values.r <- rescale(Fill.values, from = range(Fill.values), na.rm = TRUE)
  
  scale_fill_gradient2("",
                       low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 1,
                       labels = as.character(round(Breaks.fill, 3)),
                       breaks = log2(Breaks.fill), 
                       limits = range(Fill.values))
}

## -----------------------------------------------------------------------------
#' @param draws matrix of posterior draws (its x observations)
#' @param null_value numeric or vector specifying the middle value. (Defaults to 1)
# If vector it must be length(null_value) == observations
#' @sig_level numeric specifying the significance level for DPP. (Defaults to 0.6)
#' @returns list of
#				EP: probability of posterior being above null value
# 				DPP: difference in posterior probability
# 				DPP_sig: factor where 1 denotes a significant DPP
jf$getDPP <- function(draws, null_value = 1, sig_level = 0.6){
  
  if(length(null_value) != 1){
    if(length(null_value) != dim(draws)[2]){
      stop("If null_value is a vector it must be the same as the number of columns in draws. ")
    }
  }

  EP <- rep(NA, ncol(draws))
  for(i in 1:ncol(draws)){
    if(length(null_value) == 1){
      EP[i] <- mean(ifelse(draws[,i] > null_value, 1, 0))
    }else{
      EP[i] <- mean(ifelse(draws[,i] > null_value[i], 1, 0))
    }
  }
  
  # EP
  EP_low <- 1-EP
  # DPP
  DPP <- unname(abs(EP - EP_low))
  # get significant DPP
  DPP_sig = as.factor(ifelse(DPP > sig_level, 1, 0))
  # get compared_to_null
  compared_to_null = ifelse(EP > 0.8, 1, 
                           ifelse(EP < 0.2, -1, 0))
  
  # return list
  return(list(EP = EP, 
              compared_to_null = compared_to_null,
              DPP = DPP,
              DPP_sig = DPP_sig))
}

## -----------------------------------------------------------------------------
#' @param x vector of SIRs used as fill in geom_sf
#' @param extreme numeric specifying the largest SIR for which data is collapsed
#' @param nu_cats integer for the total number of labelled categories (must be even)
jf$addDiscreteColorScale <- function(x, extreme = 2, nu_cats = 6){

  nu_cats_below = nu_cats/2
  Breaks.fill <- c(min(x)-1, 
                   seq(1/extreme, 1, length = nu_cats_below)[-nu_cats_below],
                   1, 
                   seq(1, extreme, length = nu_cats_below)[-1],
                   max(x)+1)
  
  out <- cut(x, breaks = Breaks.fill)
  levels(out)[c(1,nu_cats)] <- c(paste0("<", Breaks.fill[2]),
                                  paste0(">", Breaks.fill[nu_cats]))
  
  return(out)

}

## ----------------------------------------------------------------------------
#' @param W binary weight matrx
#' @param T number of time points (defaults to NULL)
#' @returns returns a list with the following components:
# NOTE: This list is when T is defined. 
# By default `getSTcomponents(.)` returns only the BYM2 components
#		   T: total number of time points
#		   T_adj: adj object from tm.adjacency() function
#		   T_num: num object from tm.adjacency() function
#		   T_weights: num object from tm.adjacency() function
#		   L_t: length of T_adj
#		   M: total number of areas
#		   adj: adj object from as.carAdjacency() function
#		   num: num object from as.carAdjacency() function
#		   weights: weights object from as.carAdjacency() function
#		   L_s: length of adj
#		   kappa: numeric scalar for the BYM2 prior
jf$getSTcomponents <- function(W, T = NULL){
  
# get the scale and bugs components
	kappa <- jf$getBYM2scale(W)
	M_bugs <- nimble::as.carAdjacency(W)
	
if(is.null(T)){

# return the list
list(  # spatial components
	   M = nrow(W),
	   adj = M_bugs$adj,
	   num = M_bugs$num,
	   weights = M_bugs$weights,
	   L_s = length(M_bugs$adj),
	   kappa = kappa)
}else{

# Get the temporal components
T_bugs <- jf$tm.adjacency(T = T, 1)

# return the list
list(  # temporal components
	   T = T,
	   T_adj = T_bugs$T_adj,
	   T_num = T_bugs$T_num,
	   T_weights = T_bugs$T_weights,
	   L_t = length(T_bugs$T_adj),
	   # spatial components
	   M = nrow(W),
	   adj = M_bugs$adj,
	   num = M_bugs$num,
	   weights = M_bugs$weights,
	   L_s = length(M_bugs$adj),
	   kappa = kappa)
	   
}

}


## -----------------------------------------------------------------------------
#' @describe WARNING: this is a very specific function! 
jf$returnSingleYears <- function(x, fix_negs = T, OAG = F, method	= "ord"){

	rmnegs_retainsum <- function(x){
		s <- sum(x)
		x_new <- ifelse(x < 0, 0, x)
		return(s * (x_new/sum(x_new)))
	}

	temp <- graduate_beers(x$N, Age = x$age, OAG = OAG, method = method)
	if(fix_negs) temp <- rmnegs_retainsum(temp)
	data.frame(Year = x$Year[1],
			 geography_no = as.integer(x[1, 2]),
			 sex = x$sex[1],
			 age = 1:length(temp)-1,
			 N = temp)
}

## -----------------------------------------------------------------------------
#' @description Will add 8 new columns to single-age year data that groups 
# according to that specified by WA DOH
# WARNING: Very specific function (only use if `age` is a variable) 
jf$addAgeCats <- function(.data){
  .data %>% 
    mutate(age18 = cut_width(age, width = 5, boundary = 0, closed = "left", labels = F),
           age18_w = cut_width(age, width = 5, boundary = 0, closed = "left"),
           CIage4 = cut_interval(age, n = 4, breaks = c(-0.5, 17.5, 39.5, 64.5, max(age)), label = F),
           CIage4_w = cut(age, breaks = c(-0.5, 17.5, 39.5, 64.5, max(age))),
           Sage3 = cut_interval(age, n = 3, breaks = c(15.5, 44.5, 64.5, max(age)), label = F),
           Sage3_w = cut(age, breaks = c(15.5, 44.5, 64.5, max(age))),
           BoDage7 = cut_interval(age, n = 7, breaks = c(-0.5, 4.5, 14.5, 24.5, 44.5, 64.5, 84.5, max(age)), label = F),
           BoDage7_w = cut(age, breaks = c(-0.5, 4.5, 14.5, 24.5, 44.5, 64.5, 84.5, max(age))))
}

## ----------------------------------------------------------------------------
#' @description collapses an agegroup variable (which is an index) to groups of 
# size `cats_to_group` 
jf$collapseAges <- function(x, cats_to_group = 2){
  uu <- unique(x)
  mAx = uu[length(uu)]
  mIn = uu[1]
  if((length(uu)/cats_to_group) != as.integer(length(uu)/cats_to_group)){
    stop("Number of groups should be divisible!")
  }
  
  cut(x, breaks = round(length(uu)/cats_to_group), labels = FALSE)
}

## -----------------------------------------------------------------------------
#' @param adj_mat fully connected binary weight matrix
#' @return numeric scaling factor for BYM2 model
jf$getBYM2scale <- function(adj_mat){
  
  # get necessary arguments for function
  adj_mat <- as(adj_mat, "sparseMatrix")
  M <- nrow(adj_mat)
  
  # ICAR precision matrix
  Q <- Matrix::Diagonal(M, Matrix::rowSums(adj_mat)) - adj_mat
  # Add a small jitter for numerical stability
  Q_pert <- Q + Matrix::Diagonal(M) * max(Matrix::diag(Q)) * sqrt(.Machine$double.eps)
  
  # Compute the diagonal elements of the covariance matrix subject to the
  # constraint that the entries of the ICAR sum to zero
  .Q_inv <- function(Q){
    Sigma <- Matrix::solve(Q)
    A <- matrix(1,1, NROW(Sigma))
    W <- Sigma %*% t(A)
    Sigma <- Sigma - W %*% solve(A %*% W) %*% Matrix::t(W)
    return(Sigma)
  }
  
  # get Q_inv
  Q_inv <- .Q_inv(Q_pert)
  
  # Compute the geometric mean of the variances (diagonal of Q_inv)
  exp(mean(log(Matrix::diag(Q_inv))))
  
}

## -----------------------------------------------------------------------------
#' @param design_mat design matrix creates using model.matrix() function
#' @return a list with three elements:
# 		QR: list with three elements (Q_ast, R_ast, R_ast_inverse)
# 		x_c: centered (NOT scaled) design_mat with NO intercept
# 		p: number of columns of x_c
#' @details The `...` can be used to add a numeric vector of centering values that
# are passed to the center argument
jf$getQRDecomp <- function(design_mat, ...){
  # get n
  n <- nrow(design_mat)
  # center all columns
  x <- scale(design_mat, scale = F, ...)
  # remove intercept
  x <- x[,-1]
  p <- ncol(x)
  # take the QR decomposition
  qrstr <- qr(x)
  Q <- qr.Q(qrstr)
  R <- qr.R(qrstr)
  R_ast <- (1/sqrt(n - 1)) * R
  R_ast_inverse <- solve(R_ast)
  Q_ast <- Q * sqrt(n - 1)
  return(list(QR = list(Q_ast = Q_ast,
                        R_ast = R_ast,
                        R_ast_inverse = R_ast_inverse),
              x_c = x,
              p = p))
  
}

## -----------------------------------------------------------------------------
#' @param T number of time points
#' @param RW.order order of RW. Can be either 1 or 2
#' @returns a list containing three components: adj, num and weights
# Functions stolen from https://www.sptmbook.com/datacode.html
jf$tm.adjacency <- function(T,RW.order) {
  if (RW.order==1) {
    ####   for random walk order 1
    num <- numeric(T)
    weights <- numeric((T-2)*2+2)
    adj <- numeric((T-2)*2+2)
    for (t in 1:1) {
      weights[t] <- 1
      adj[t] <- t+1
      num[t] <- 1
    }
    for (t in 2:(T-1)) {
      weights[2+(t-2)*2] <- 1
      adj[2+(t-2)*2] <- t-1
      weights[3+(t-2)*2] <- 1
      adj[3+(t-2)*2] <- t+1
      num[t] <- 2
    }
    for (t in T:T) {
      weights[(T-2)*2+2] <- 1
      adj[(T-2)*2+2] <- t-1
      num[T] <- 1
    }
  }
  if (RW.order==2) {
    ####   for random walk order 2
    num <- numeric(T)
    weights <- numeric((T-4)*4+3*2+2*2)
    adj <- numeric((T-4)*4+3*2+2*2)
    for (t in 1:1) {
      weights[t] <- 2
      adj[t] <- t+1
      weights[t+1] <- -1
      adj[t+1] <- t+2
      num[t] <- 2
    }
    for (t in 2:2) {
      weights[t+1] <- 2
      adj[t+1] <- t-1
      weights[t+2] <- 4
      adj[t+2] <- t+1
      weights[t+3] <- -1
      adj[t+3] <- t+2
      num[t] <- 3
    }
    for (t in 3:(T-2)) {
      weights[6+(t-3)*4] <- -1
      adj[6+(t-3)*4] <- t-2
      weights[7+(t-3)*4] <- 4
      adj[7+(t-3)*4] <- t-1
      weights[8+(t-3)*4] <- 4
      adj[8+(t-3)*4] <- t+1
      weights[9+(t-3)*4] <- -1
      adj[9+(t-3)*4] <- t+2
      num[t] <- 4
    }
    for (t in (T-1):(T-1)) {
      weights[(T-4)*4+6] <- 2
      adj[(T-4)*4+6] <- t+1
      weights[(T-4)*4+7] <- 4
      adj[(T-4)*4+7] <- t-1
      weights[(T-4)*4+8] <- -1
      adj[(T-4)*4+8] <- t-2
      num[t] <- 3
    }
    for (t in T:T) {
      weights[(T-4)*4+9] <- 2
      adj[(T-4)*4+9] <- t-1
      weights[(T-4)*4+10] <- -1
      adj[(T-4)*4+10] <- t-2
      num[t] <- 2
    }
  }  
  if (T==2) {
    adj <- c(2,1)
    weights <- c(1,1)
    num <- c(1,1)
  }
  adj.tm <- list()
  adj.tm$adj <- adj
  adj.tm$num <- num
  adj.tm$weights <- weights
  
  # Add index for names
  names(adj.tm) <- paste0("T_", names(adj.tm))
  return(adj.tm)
}

## -----------------------------------------------------------------------------
#' @param y vector of non-integer counts
#' @param N vector of corresponding populations
#' @description Function used to calculate rounded standardised counts for
# the area model for ASRs
#' @return Returns a dataframe of `length(y)` with (y, N, y_tilde, N_tilde)
jf$sIntRound <- function(y, N){

data.frame(y = y) %>%
	mutate(# rounded integers using ceiling
		   y_tilde = ceiling(y),
		   # adjustment factor
		   cc = y_tilde/y,
		   cc = ifelse(y == 0,1,cc),
		   # adjusted population
		   N_tilde = ifelse(y == 0, N, N * cc),
		   # original population
		   N = N)
}

## -----------------------------------------------------------------------------
#' @param x numeric vector that should be scaled. Most likely `y_tilde`
#' @param cc column given as part of the output of `jf$sIntRound`
#' @returns a numeric vector of the same length as x
#' @details Can be used in `apply` or `pbapply` code. 
jf$sIntRound_reverse <- function(x, cc){
  stopifnot(length(x) == length(cc))
  x/cc
}


## -----------------------------------------------------------------------------
#' @description requires the following packages; `nimble` and `posterior`
#' @param code
#' @param nD data in list format
#' @param nC constants in list format
#' @param nI function that generates initial values
#' @param monitors character vector of parameters to sample
#' @param niter number of samples
#' @param nburnin number of sample to be discarded (defaults to niter/2)
#' @param thin frequency of thinning (defaults to 1)
#' @param nchains (defaults to 4)
#' @param onlySlice logical for the mcmcConfig() function (defaults to F)
#' @param print_samplers logical for whether to display the samplers used (defaults to F)
#' @param optimBeta Logical to set the sampler for regression coefficients (defaults to F)
#' @param beta_name character for the variable name of the regression coefficients (default to NULL)
#' @param sampler_name character defining the sampler for regression coefficients (defaults "RW_block")
#' @param adaptInterval frequency of adaption in burnin (defaults to 200)
#' @param Cmcmc previous Cmcmc object (useful when a model should be rerun with more iterations)
# (defaults to NULL)
#' @param verbose prints convergence diagnostics after run (defaults to F)
#' @return a list of four elements: 
#				fit: original Nimble fit object
# 				specs: list of MCMC settings such as runtime, niter,nburnin,thin,nchains
#				summary: summary object of monitored parameters (includes Rhat and ESS)
#				message: copy of the convergence messages printed on completion. 
#						 To redisplay the convergence messages type `message(fit$messages)` 
#				MCMCMsettings: data.frame displaying the MCMCsettings including runTime (mins, hours, days)
# 							   and average posterior draws per minute and the number of useable draws
# 				Cmcmc: The compiled nimble object (useful for reruns of the same model 
# 					   where recompiling is not necessary)
#' @description Number of sampler kept per chain is (niter-nburnin)/thin
# To define a `ess` univariate sampler for a vector, set `beta_name` to `paste0("B_qr[", 1:16, "]")` for example
jf$SampleNimble <-function(code, 
                           nD, 
                           nC, 
                           nI, 
                           monitors, 
                           niter,
                           nburnin = niter/2, 
                           thin = 1, 
                           nchains = 4,
						   print_samplers = FALSE,
						   onlySlice = FALSE,
						   adaptInterval = 200,
                           optimBeta = FALSE,
                           beta_name = NULL,
                           sampler_name = "RW_block",
						   Cmcmc = NULL,
						   verbose = FALSE,
						   seed = 45){
						   
	# If using previous Cmcmc object then immediately start sampling
	if(is.null(Cmcmc)){
  
	# data lists
	Rmodel <- nimbleModel(code = code, 
						data = nD,
						inits = nI(),
						constants = nC)
	Rmodel$initializeInfo()

	# generate a default configuration for the MCMC
	mcmcConf <- configureMCMC(Rmodel, 
							monitors = monitors,
							enableWAIC = T, 
							print = F,
							control = list(adaptInterval = adaptInterval),
							onlySlice = onlySlice)
	# use RW for betas
	if(optimBeta == T){
		if(is.null(beta_name)){
		  stop("Please provide the name of the beta vector!")
		}
		if(length(beta_name) > 1){
		  for(i in 1:length(beta_name)){
			mcmcConf$removeSampler(beta_name[i])
			mcmcConf$addSampler(type = sampler_name, target = beta_name[i])
		  }
		}else{
		  mcmcConf$removeSampler(beta_name)
		  mcmcConf$addSampler(type = sampler_name, target = beta_name)
		}
		# print samplers
		print(mcmcConf) 
	}else{
		if(print_samplers){
			print(mcmcConf)
			# Interim check
			  if(readline(prompt = "To continue please type 'Y' > ") != "Y"){
				stop("Stopping...")
			  }
		}
	}

	# build and compile the MCMC
	Rmcmc <- buildMCMC(mcmcConf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
	
	} # If Cmcmc object supplied then start function here
  
	# sample
	message("Starting sampling for ", (niter - nburnin)/thin, " iterations for each of ", nchains, " chains.")
	m_s <- Sys.time()
	set.seed(seed)	
	fit <- runMCMC(Cmcmc, 
				 inits = nI(),
				 niter = niter,
				 nburnin = nburnin, 
				 thin = thin, 
				 nchains = nchains, 
				 summary = T,
				 WAIC = T,
				 samplesAsCodaMCMC = T)
	time <- as.numeric(Sys.time() - m_s, units = "mins")
	message("Sampling took ", round(time, 2), " mins")
	
	# Resolve issue with columns of summ being of class
	# "pillar_num"  "pillar_vctr" "vctrs_vctr"  "double"  
	is.pillar_num <- function(x){
	  class(x)[1] == "pillar_num"
	}

	# create summary and get metrics
	summ <- fit$samples %>% 
	posterior::summarise_draws() %>% 
	mutate(eb_per_sec = ess_bulk/(time*60),
	       et_per_sec = ess_tail/(time*60)) %>% 
	  mutate_if(is.pillar_num, as.numeric)
  
	# report convergence messages
	minimum_ess <- nchains * 100
	mcmc_diags <- paste0("Median Rhat: ", round(median(summ$rhat, na.rm = T), 3), " \n",
				round(100*mean(ifelse(summ$rhat > 1.01, 1, 0), na.rm = T), 2), "% of Rhats larger than 1.01 \n",
				"Max Rhat = ", round(max(summ$rhat, na.rm = T), 2), " (", summ$variable[which(summ$rhat == max(summ$rhat, na.rm = T))][1], ") \n",
				round(100*mean(ifelse(summ$ess_bulk < minimum_ess, 1, 0), na.rm = T), 2), "% of ess_bulk are too small \n", 
				"Min ess_bulk = ", round(min(summ$ess_bulk, na.rm = T), 2), " (", summ$variable[which(summ$ess_bulk == min(summ$ess_bulk, na.rm = T))][1], ") \n",
				round(100*mean(ifelse(summ$ess_tail < minimum_ess, 1, 0), na.rm = T), 2), "% of ess_tail are too small \n",
				"Min ess_tail = ", round(min(summ$ess_tail, na.rm = T), 2), " (", summ$variable[which(summ$ess_tail == min(summ$ess_tail, na.rm = T))][1], ") \n",
				"Average posterior draws per minute: ", round((niter*nchains)/time, 2))

	# report messages
	if(verbose){
		message(mcmc_diags)
	}

	# MCMC settings
	MCMCsettings <- data.frame(metric = c("niter", "nburnin", "thin", "nchains", 
										"Useable draws", "av_pd_m", 
										"runTime_mins", "runTime_hours", "runTime_days"),
								value = c(niter, nburnin, thin, nchains, 
										 nchains*((niter-nburnin)/thin), round((niter*nchains)/time, 2),
										 round(time, 2), round(time/60, 2), round(time/(60*24), 2)))

	# return objects
	return(list(fit = fit,
				specs = list(time = time, 
				niter = niter,
				nburnin = nburnin, 
				thin = thin, 
				nchains = nchains),
				summary = summ,
				messages = mcmc_diags,
				MCMCsettings = MCMCsettings,
				Cmcmc = Cmcmc))
}

## -----------------------------------------------------------------------------
#' @description requires the following packages; `nimble` and `posterior`
#' @param code
#' @param nD data in list format
#' @param nC constants in list format
#' @param nI function that generates initial values
#' @param monitors character vector of parameters to sample
#' @param niter number of samples
#' @param nburnin number of sample to be discarded (defaults to niter/2)
#' @param thin frequency of thinning (defaults to 1)
#' @param nchains (defaults to 4)
#' @param onlySlice logical for the mcmcConfig() function (defaults to F)
#' @param print_samplers logical for whether to display the samplers used (defaults to F)
#' @param optimBeta Logical to set the sampler for regression coefficients (defaults to F)
#' @param beta_name character for the variable name of the regression coefficients (default to NULL)
#' @param sampler_name character defining the sampler for regression coefficients (defaults "RW_block")
#' @param adaptInterval frequency of adaption in burnin (defaults to 200)
#' @param Cmcmc previous Cmcmc object (useful when a model should be rerun with more iterations)
# (defaults to NULL)
#' @return a list of four elements: 
#				fit: original Nimble fit object
#				specs: list of MCMC settings such as runtime, niter,nburnin,thin,nchains
# 			Cmcmc: The compiled nimble object (useful for reruns of the same model 
# 					   where recompiling is not necessary)
#' @description Number of sampler kept per chain is (niter-nburnin)/thin
# To define a `ess` univariate sampler for a vector, set `beta_name` to `paste0("B_qr[", 1:16, "]")` for example
jf$runNimble <-function(code, 
                           nD, 
                           nC, 
                           nI, 
                           monitors, 
                           niter,
                           nburnin = niter/2, 
                           thin = 1, 
                           nchains = 4,
						   print_samplers = FALSE,
						   onlySlice = FALSE,
						   adaptInterval = 200,
                           optimBeta = FALSE,
                           beta_name = NULL,
                           sampler_name = "RW_block",
						   Cmcmc = NULL,
						   verbose = TRUE,
						   seed = 45){
						   
	# If using previous Cmcmc object then immediately start sampling
	if(is.null(Cmcmc)){
  
	# data lists
	Rmodel <- nimbleModel(code = code, 
						data = nD,
						inits = nI(),
						constants = nC)
	Rmodel$initializeInfo()

	# generate a default configuration for the MCMC
	mcmcConf <- configureMCMC(Rmodel, 
							monitors = monitors,
							enableWAIC = T, 
							print = F,
							control = list(adaptInterval = adaptInterval),
							onlySlice = onlySlice)
	# use RW for betas
	if(optimBeta == T){
		if(is.null(beta_name)){
		  stop("Please provide the name of the beta vector!")
		}
		if(length(beta_name) > 1){
		  for(i in 1:length(beta_name)){
			mcmcConf$removeSampler(beta_name[i])
			mcmcConf$addSampler(type = sampler_name, target = beta_name[i])
		  }
		}else{
		  mcmcConf$removeSampler(beta_name)
		  mcmcConf$addSampler(type = sampler_name, target = beta_name)
		}
		# print samplers
		print(mcmcConf) 
	}else{
		if(print_samplers){
			print(mcmcConf)
			# Interim check
			  if(readline(prompt = "To continue please type 'Y' > ") != "Y"){
				stop("Stopping...")
			  }
		}
	}

	# build and compile the MCMC
	Rmcmc <- buildMCMC(mcmcConf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
	
	} # If Cmcmc object supplied then start function here
  
	# sample
	message("Starting sampling for ", (niter - nburnin)/thin, " iterations for each of ", nchains, " chains.")
	m_s <- Sys.time() 
	set.seed(seed)
	fit <- runMCMC(Cmcmc, 
				 inits = nI(),
				 niter = niter,
				 nburnin = nburnin, 
				 thin = thin, 
				 nchains = nchains, 
				 summary = T,
				 WAIC = T,
				 samplesAsCodaMCMC = T)
	time <- as.numeric(Sys.time() - m_s, units = "mins")
	message("Sampling took ", round(time, 2), " mins")

	# return objects
	return(list(fit = fit, 
              specs = list(time = time, 
                           niter = niter,
                           nburnin = nburnin, 
                           thin = thin, 
                           nchains = nchains),
				      Cmcmc = Cmcmc))
}

## -----------------------------------------------------------------------------
#' @description requires the following packages; `nimble` and `posterior`
#' @param input input list from `runNimble`
#' @param trim_its (defaults to NULL) numeric defining how many iterations to discard
#' @param verbose prints convergence diagnostics after run (defaults to T)
#' @return a list of four elements: 
#				fit: original Nimble fit object
#				summary: summary object of monitored parameters (includes Rhat and ESS)
#				message: copy of the convergence messages printed on completion. 
#						 To redisplay the convergence messages type `message(fit$messages)` 
#				MCMCMsettings: data.frame displaying the MCMCsettings including runTime (mins, hours, days)
# 							   and average posterior draws per minute and the number of useable draws
jf$sumNimble <-function(input, trim_its = NULL, verbose = TRUE){
						   
  time = input$specs$time
  niter = input$specs$niter
  nburnin = input$specs$nburnin
  thin = input$specs$thin
  nchains = input$specs$nchains
  ud <- nchains*((niter-nburnin)/thin)
	
	# Resolve issue with columns of summ being of class
	# "pillar_num"  "pillar_vctr" "vctrs_vctr"  "double"  
	is.pillar_num <- function(x){
	  class(x)[1] == "pillar_num"
	}
	
	if(!is.null(trim_its)){
		temporary_samples <- input$fit$samples
		for(j in 1:length(temporary_samples)){
			temporary_samples[[j]] <- temporary_samples[[j]][-c(1:trim_its),]
		}
		input$fit$samples <- temporary_samples
		message(paste0("Removing the first ", trim_its, " draws..."))

		# update iterations
		nburnin = nburnin + trim_its
		ud = ud - (nchains * trim_its)
	}

	# create summary and get metrics
	summ <- input$fit$samples %>% 
	posterior::summarise_draws() %>% 
	mutate(eb_per_sec = ess_bulk/(time*60),
	       et_per_sec = ess_tail/(time*60)) %>% 
	  mutate_if(is.pillar_num, as.numeric)
  
	# report convergence messages
	minimum_ess <- input$specs$nchains * 100
	mcmc_diags <- paste0("Median Rhat: ", round(median(summ$rhat, na.rm = T), 3), " \n",
				round(100*mean(ifelse(summ$rhat > 1.01, 1, 0), na.rm = T), 2), "% of Rhats larger than 1.01 \n",
				"Max Rhat = ", round(max(summ$rhat, na.rm = T), 2), " (", summ$variable[which(summ$rhat == max(summ$rhat, na.rm = T))][1], ") \n",
				round(100*mean(ifelse(summ$ess_bulk < minimum_ess, 1, 0), na.rm = T), 2), "% of ess_bulk are too small \n", 
				"Min ess_bulk = ", round(min(summ$ess_bulk, na.rm = T), 2), " (", summ$variable[which(summ$ess_bulk == min(summ$ess_bulk, na.rm = T))][1], ") \n",
				round(100*mean(ifelse(summ$ess_tail < minimum_ess, 1, 0), na.rm = T), 2), "% of ess_tail are too small \n",
				"Min ess_tail = ", round(min(summ$ess_tail, na.rm = T), 2), " (", summ$variable[which(summ$ess_tail == min(summ$ess_tail, na.rm = T))][1], ") \n",
				"Average posterior draws per minute: ", round((niter*nchains)/time, 2))

	# report messages
	if(verbose){
		message(mcmc_diags)
	}

	# MCMC settings
	MCMCsettings <- data.frame(metric = c("niter", "nburnin", "thin", "nchains", 
										"Useable draws", "av_pd_m", 
										"runTime_mins", "runTime_hours", "runTime_days"),
								value = c(niter, nburnin, thin, nchains, 
										 ud, round((niter*nchains)/time, 2),
										 round(time, 2), round(time/60, 2), round(time/(60*24), 2)))

	# return objects
	return(list(fit = input$fit, 
				summary = summ,
				messages = mcmc_diags,
				MCMCsettings = MCMCsettings))
}

## -----------------------------------------------------------------------------
#' @param draws iterations by observations
#' @param model model column content (as single character)
#' @param metric metric column content (as single character)
#' @param conf_level numeric for the confidence size of the interval (defaults to 0.95)
#' @param other_data data.frame with other data to be column binded to result 
# (NOTE: must be same dimensions and same order!)
#' @param prefix character vector added to point, upper, lower, se and RSE columns (defaults to "")
#' @param addDPP logical (defaults to false)
#' @returns Dataset with posterior summaries including: 
# median point estimates, credible intervals (HDI)
# standard deviations and RSE
# More arguments can be passed to the getDPP() function using `...` (e.g. null_value)
jf$getResultsData <- function(draws, 
                           model = NULL, metric = NULL,
                           prefix = "",
                           conf_level = 0.95,
                           other_data = NULL, addDPP = FALSE, ...){
  
  # sum_func <- function(x){
  #   c(point = median(x, na.rm = T),
  #     lower = unname(HDInterval::hdi(x, credMass = conf_level)[1]),
  #     upper = unname(HDInterval::hdi(x, credMass = conf_level)[2]),
  #     se = sd(x, na.rm = T))
  # }
  # bind_rows(pblapply(asplit(draws, 2), sum_func))
  
  if(!is.null(other_data)){
    message(paste0("NOTE (not an error): Please check that the row order of ", deparse(substitute(other_data)), " matches that of the column order of ", deparse(substitute(draws))))
    if(nrow(other_data) != ncol(draws))stop(paste0("The number of columns of ", deparse(substitute(draws)), " does NOT match the number of rows of ", deparse(substitute(other_data)), ". They must match!"))
  }
  
  if(is.null(dim(draws))){
    r <- data.frame(point = draws,
                    lower = NA,
                    upper = NA,
                    se = NA) %>%
      mutate(RSE = NA) %>%
      setNames(paste0(prefix, names(.)))
    r <- bind_cols(r, other_data)
  }else{
    # Get objects
    message("Progress ... -> Point estimates...")
    point_in = pbapply::pbapply(draws, 2, median, na.rm = T)
    message("Progress ... -> Standard errors...")
    se_in = pbapply::pbapply(draws, 2, sd, na.rm = T)
    message("Progress ... -> Highest density intervals...")
    hd_ints = pbapply::pbapply(draws, 2, HDInterval::hdi, credMass = conf_level)
    if(addDPP){
      DPP <- jf$getDPP(draws, ...)
      r <- data.frame(point = point_in,
                      #lower = apply(draws, 2, quantile, prob = 0.025, na.rm = T),
                      lower = hd_ints[1,],
                      #upper = apply(draws, 2, quantile, prob = 0.975, na.rm = T),
                      upper = hd_ints[2,],
                      se = se_in,
                      EP = DPP$EP,
                      compared_to_null = DPP$compared_to_null,
                      DPP = DPP$DPP,
                      DPPsig = DPP$DPP_sig) %>%
        mutate(RSE = 100 * (se/point)) %>%
        setNames(paste0(prefix, names(.)))
      r <- bind_cols(r, other_data)
    }else{
      r <- data.frame(point = point_in,
                      #lower = apply(draws, 2, quantile, prob = 0.025, na.rm = T),
                      lower = hd_ints[1,],
                      #upper = apply(draws, 2, quantile, prob = 0.975, na.rm = T),
                      upper = hd_ints[2,],
                      se = se_in) %>%
        mutate(RSE = 100 * (se/point)) %>%
        setNames(paste0(prefix, names(.)))
      r <- bind_cols(r, other_data)
    }
  }
  
  # Add columns if given
  if(!is.null(model) & !is.null(metric)){
    r <- r %>% 
      mutate(model = model,
             metric = metric) %>% 
      relocate(model, metric)
  }
  if(!is.null(model) & is.null(metric)){
    r <- r %>% 
      mutate(model = model) %>% 
      relocate(model)
  }
  if(is.null(model) & !is.null(metric)){
    r <- r %>% 
      mutate(metric = metric) %>% 
      relocate(metric)
  }
  
  # return objects
  rownames(r) <- NULL
  return(r)
}

## -----------------------------------------------------------------------------
#' @param draws any matrix of posterior draws (iterations must be rows)
#' @param thin numeric amount to thin
jf$postThin <- function(draws, thin){
  thin_seq <- seq(1, nrow(draws), thin)
  thinned_draws <- draws[thin_seq,]
  message(paste0(nrow(draws), " draws given. ", nrow(thinned_draws), " draws returned."))
  return(thinned_draws)
}

## -----------------------------------------------------------------------------
#' @param y vector of binary outcomes (length N)
#' @param dom vector of indexes for areas (length N)
#' @param sweight vector of sample weights (length N)
#' @param domsize dataframe (M x 2): first column - area indexes, 2nd column: population
#' @param area_name character for the variable that specifies the subgroups (defaults to "area")
#' @param eps perturbation amount (defaults to 0.001 for estimates of zero)
#' @param verbose logical (defaults to F). Reserved for debugging
#' @returns dataframe of length M with a variety of area-level metrics
jf$jDirect <- function(y, dom, sweight, domsize, area_name = "area", eps = 0.001, verbose = F){
  
  result <- data.frame(Domain=0,SampSize=0,Direct=0,
                       VAR=0,SD=0,CV=0)
  
  did     <- unique(dom)    # unique identifiers of domains
  Dsample <- length(did)    # number of domains in sample
  
  # Calculate HT direct estimator for sampled domains   
  nds      <-rep(0,Dsample)   # domain sample sizes
  dirds    <-rep(0,Dsample)   # domain direct estimators 
  vardirds <-rep(0,Dsample)   # variances of direct estimators
  N        <-rep(0,Dsample)   # population size
  
  for (d in 1:Dsample){
    if(verbose){
      message("Starting the ", d, "th area.")
    }
    yd       <- y[dom==did[d]]
    nds[d]   <- length(yd)
    
    sweightd <- sweight[dom==did[d]]
    w <- nds[d] * ( sweightd / sum(sweightd) )
    N[d] <- unlist(domsize[(domsize[,1]==did[d]),2])
    
    dirds[d] <- (sum(yd*w)/sum(w))[[1]] 
    
    # Approximated unbiased estimator of variance of HT direct estimator
    vardirds[d]<-(1/nds[d]) * (1-(nds[d]/N[d])) * (1/(nds[d]-1)) * sum( w^2 * (yd - dirds[d])^2 )
  }
  
  # create output dataset
  out <- data.frame(area = did, 
                    n = nds, 
                    N = N,
                    phat = dirds,
                    phat_VAR = vardirds) %>% 
    mutate(# create unstable indicator
           unstable = phat %in% c(0,1),
           # variance quantities
           phat_VAR = ifelse(unstable, NA, phat_VAR),
		   phat_SE = sqrt(phat_VAR),
		   phat_RSE = 100 * (phat_SE/phat),
		   # add intervals for phat estimates
		   phat_lower = ifelse(phat - 1.96 * phat_SE < 0, 0, 
							ifelse(phat - 1.96 * phat_SE > 1, 1, phat - 1.96 * phat_SE)),
		   phat_upper = ifelse(phat + 1.96 * phat_SE < 0, 0, 
							ifelse(phat + 1.96 * phat_SE > 1, 1, phat + 1.96 * phat_SE)),
           # add small perturbation to phat
           phat = ifelse(phat == 0, eps, phat),
           phat = ifelse(phat == 1, 1-eps, phat),
           # create empirical logit
           phat_u = jf$jlogit(phat),
           phat_u_VAR = phat_VAR/((phat * (1 - phat))^2),
           phat_u_SE = sqrt(phat_u_VAR),
		   # create y_tilde and n_tilde
		   n_tilde = (phat * (1 - phat))/phat_VAR,
		   y_tilde = n_tilde * phat)
	names(out)[1] <- area_name
	
  # return the dataframe
  return(out)
}

## -----------------------------------------------------------------------------
# Depreciated function for CAR models in Stan
jf$mungeCARdata4stan = function(adjBUGS,numBUGS) {
  N = length(numBUGS);
  nn = numBUGS;
  N_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=N_edges);
  node2 = vector(mode="numeric", length=N_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:N) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("N" = N,
			   "N_edges" = N_edges,
				"node1" = node1, 
				"node2" = node2));
}

## ----------------------------------------------------------------------------
#' @param x1 vector of posterior draws
#' @param x2 vector of the same length as x1
#' @param x3 vector of the same length as x1
#' @return density plot comparing x1 and x2
jf$compareDensity <- function(x1, x2, x3 = NULL){
	if(is.null(x3)){
	data.frame(x1 = x1, x2 = x2) %>% 
	pivot_longer(everything()) %>% 
		ggplot(aes(x = value, col = name, fill = name))+
			geom_density(alpha = 0.2)+theme_bw()+
			labs(subtitle = paste0("Median (sd): x1 - ", 
								 round(median(x1),2), "(", round(sd(x1),2), ") ",
								 "x2 - ", round(median(x2),2), " (", round(sd(x2),2),")"))
	}else{
	data.frame(x1 = x1, x2 = x2, x3 = x3) %>% 
	pivot_longer(everything()) %>% 
		ggplot(aes(x = value, col = name, fill = name))+
			geom_density(alpha = 0.2)+theme_bw()+
			labs(subtitle = paste0("Median (sd): x1 - ", 
								 round(median(x1),2), "(", round(sd(x1),2), ") ",
								 "x2 - ", round(median(x2),2), " (", round(sd(x2),2),")",
								 "x3 - ", round(median(x3),2), " (", round(sd(x3),2),")"))
	}
}

## ----------------------------------------------------------------------------
#' @param av_pd_m average number of draws per minute
#' @param available_hours
#' @param ud numeric vector specifying the total number of useful draws we need
#' @param nchains defaults to 4
#' @param thin defaults to 1
#' @param defaults to half of niter
#' @return by entering `av_pd_m` and `avail_hours`, the function returns the best MCMC settings
#' @return by entering `av_pd_m` and `ud`, the function returns an approximation of the run time
jf$MCMCrecommendations <- function(av_pd_m,
									avail_hours = NULL,
									ud = NULL,
									nchains = 4, 
									thin = 1,
									nburnin = NULL){
  
  # estimate draws based on speed and available time
  if(is.null(ud)){
    total_mins <- avail_hours * 60
    total_its <- total_mins*av_pd_m
    niter <- total_its/nchains
    if(is.null(nburnin)){
      nburnin <- niter/2
    }
    ud <- ((niter - nburnin)/thin)*nchains
    
    data.frame(metric = c("nchains", "thin", "niter", "nburnin", "Useful draws"),
               recommendations = round(c(nchains, thin, niter, nburnin, ud))) 
  }else{
  # estimate run time based on speed and niter, etc
    if(is.null(nburnin)){
      nburnin <- (ud/nchains)*thin
    }
    niter <- ((ud/nchains)*thin)+(nburnin)
	mins <- (niter*nchains)/av_pd_m
    hours <- mins/60
    days <- hours/24
    
    data.frame(metric = c("nchains", "thin", "niter", "nburnin", "Useful draws", "mins", "hours", "days"),
               recommendations = c(round(c(nchains, thin, niter, nburnin, ud)), round(c(mins, hours, days), 2)))
  }
}

## --------------------------------------------------------------------------------------
#' @description requires the following packages; `nimble` and `posterior`
#' @param formula formula for the ST.CARanova() function
#' @param burnin number of samples for burnin (defaults to half of n.sample)
#' @param n.sample total number of samples (defaults to 20,000)
#' @param thin number of samples to thin (defaults to 10)
#' @param data dataset to pass to ST.CARanova()
#' @param W binary weight matrix
#' @param area character for area variable
#' @param year character for year variable
#' @param ofs numeric vector of offset term (for example, E or N_tilde)
#' @param y numeric vector of counts
#' @param prior.tau2 two dimension vector specifying the two parameters of
# the inverse gamma prior on the variance parameters (defaults to `c(1, 0.01)`)
#' @param post_integer_adj_scaling (defaults to NULL) Numeric vector of scaling factors to be applied to fitted counts
# Transformation is applied after residuals are derived. 
# Should be the `cc` column provided by the `jf$sIntRound()` function.  
#' @param MALA logical (defaults to FALSE) for whether to use MALA for fixed effects
#' @param verbose logical (defaults to TRUE) to print the chain progress
#' @param interaction include space-time random effect (defaults to TRUE)
#' @param ... arguments to the ST.CARanova() function can be used
#' @details Automatically sets `rho.T = T` and `family = 'poisson'`
# The function will automatically order the data correctly. 
# NOTE: Will always use 4 chains!
#' @returns list with the following:
# 			fit: list where each element is a CARBayesST object
#			draws: list with mcmc.list objects for each parameter
#			sres_draws: matrix (array) of standardized residuals (iterations by obs)
#			fitted_draws: matrix (array) of draws of fitted values (iterations by obs)
#			rate_draws: matrix (array) of draws of rate values (fitted values divided by ofs) (iterations by obs)
#			summary: summary of all model parameters across chains
#			messages: use message(fit$messages) to get convergence messages
#			modelfit: modelfit metrics for each chain 			  
#			time: runtime in mins
jf$SampleCBST <- function(formula,
                          data, 
                          W,
                          area,
                          year,
                          ofs,
                          y,
                          burnin = n.sample/2, 
                          n.sample = 20000, 
                          thin = 10,
                          prior.tau2 = c(1, 0.01),
						  post_integer_adj_scaling = NULL, 
                          MALA = FALSE,
                          verbose = TRUE,
                          interaction = TRUE){
  
  # Message
  ud <- 4*( (n.sample-burnin)/thin )
  if(verbose) message("Useable draws: ",  ud)
  
  # Ensure ordering is correct
  foo <- function(.data, year, area){
    .data %>% 
      arrange(.data[[year]], .data[[area]])
  }
  temp1 <- data %>% foo(year, area)
  if(!identical(head(data), head(temp1))){
    stop("Please check the order of your input data.\nThe first K observations should be for\ntime point 1!")
  }
  
  # check dimension of scaling factor
  if(!is.null(post_integer_adj_scaling)){
	stopifnot(length(post_integer_adj_scaling) == nrow(data))
  }
  
  # Fit models
  s <- Sys.time()
  fit <- replicate(4, 
                   ST.CARanova(as.formula(formula), 
                               family = "poisson", 
                               burnin = burnin, 
                               n.sample = n.sample, 
                               thin = thin,
                               rho.T = 1,
                               data = data, 
                               W = W,
                               trials = NULL,
                               interaction = interaction,
                               MALA = MALA,
                               verbose = verbose), 
                                simplify = F)
  time <- as.numeric(Sys.time() - s, units = "mins")
  if(verbose) message("Sampling took ",  round(time, 2), " mins.")
  
  # get modelfit
  modelfit <- list(fit[[1]]$modelfit, 
                   fit[[2]]$modelfit,
                   fit[[3]]$modelfit,
                   fit[[4]]$modelfit) %>% 
    bind_rows() %>% 
    mutate(chain = 1:4) %>%
    relocate(chain)
  
  # get fitted draws
  count_draws <- as.matrix(coda::mcmc.list(fit[[1]]$samples$fitted, 
                                           fit[[2]]$samples$fitted, 
                                           fit[[3]]$samples$fitted, 
                                           fit[[4]]$samples$fitted))
										   
  # get rate draws
  rate_draws <- t(t(unname(count_draws))/ofs)
  
  # Get standardized residuals
  y_mat <- matrix(rep(y, ud/4), byrow = T, nrow = ud/4)
  res1 <- (y_mat - fit[[1]]$samples$fitted)/sqrt(fit[[1]]$samples$fitted)
  res2 <- (y_mat - fit[[2]]$samples$fitted)/sqrt(fit[[2]]$samples$fitted)
  res3 <- (y_mat - fit[[3]]$samples$fitted)/sqrt(fit[[3]]$samples$fitted)
  res4 <- (y_mat - fit[[4]]$samples$fitted)/sqrt(fit[[4]]$samples$fitted)
  sres_draws <- unname(as.matrix(coda::mcmc.list(res1, res2, res3, res4)))
  
  # apply post integer adjustment
  if(!is.null(post_integer_adj_scaling)){
	count_draws <- t(apply(count_draws, 1, jf$sIntRound_reverse, cc = post_integer_adj_scaling))
  }
  fitted_draws <- unname(count_draws)
  
  # Create summary dataframe
  # beta
  beta_cats <- ncol(fit[[1]]$samples$beta)
  for(i in 1:4){
    attr(fit[[i]]$samples$beta, "dimnames") <- list(NULL, paste0("beta[", rep(1:beta_cats), "]"))
  }
  beta_draws <- coda::mcmc.list(fit[[1]]$samples$beta, 
                                fit[[2]]$samples$beta,
                                fit[[3]]$samples$beta, 
                                fit[[4]]$samples$beta)
  beta <- as.matrix(beta_draws) %>% 
    posterior::summarise_draws() %>% 
    mutate(variable = paste0("beta[", 1:nrow(.), "]"),
           eb_per_sec = ess_bulk/(time*60),
           et_per_sec = ess_tail/(time*60))
  
  # tau -> we call them sigma
  if(interaction){char_vec <- paste0("sigma2_", c("theta", "gamma", "delta"))
  }else{char_vec <- paste0("sigma2_", c("theta", "gamma"))}
  for(i in 1:4){
    attr(fit[[i]]$samples$tau2, "dimnames")[[2]] <- char_vec
  }
  sigma2_draws <- coda::mcmc.list(fit[[1]]$samples$tau2, 
                                  fit[[2]]$samples$tau2,
                                  fit[[3]]$samples$tau2,
                                  fit[[4]]$samples$tau2)
  #sigma2_draws <- coda::mcmc.list(lapply(fit, function(x) x$samples$tau2))
  sigma2 <- as.matrix(sigma2_draws) %>% 
    posterior::summarise_draws() %>% 
    mutate(variable = char_vec,
           eb_per_sec = ess_bulk/(time*60),
           et_per_sec = ess_tail/(time*60))
  
  # rho
  rho_draws <- coda::mcmc.list(fit[[1]]$samples$rho, 
                               fit[[2]]$samples$rho,
                               fit[[3]]$samples$rho,
                               fit[[4]]$samples$rho)
  rho <- as.matrix(rho_draws) %>% 
    posterior::summarise_draws() %>% 
    mutate(variable = "rho_S",
           eb_per_sec = ess_bulk/(time*60),
           et_per_sec = ess_tail/(time*60))
  
  # phi -> spatial random effects (theta)
  theta_cats <- ncol(fit[[1]]$samples$phi)
  for(i in 1:4){
    attr(fit[[i]]$samples$phi, "dimnames") <- list(NULL, paste0("theta[", rep(1:theta_cats), "]"))
  }
  theta_draws <- coda::mcmc.list(fit[[1]]$samples$phi, 
                                 fit[[2]]$samples$phi,
                                 fit[[3]]$samples$phi,
                                 fit[[4]]$samples$phi)
  theta <- as.matrix(theta_draws) %>% 
    posterior::summarise_draws() %>% 
    mutate(variable = paste0("theta[", 1:nrow(.), "]"),
           eb_per_sec = ess_bulk/(time*60),
           et_per_sec = ess_tail/(time*60))
  
  # delta -> temporal random effects (gamma)
  gamma_cats <- ncol(fit[[1]]$samples$delta)
  for(i in 1:4){
    attr(fit[[i]]$samples$delta, "dimnames") <- list(NULL, paste0("gamma[", rep(1:gamma_cats), "]"))
  }
  gamma_draws <- coda::mcmc.list(fit[[1]]$samples$delta, 
                                 fit[[2]]$samples$delta,
                                 fit[[3]]$samples$delta,
                                 fit[[4]]$samples$delta)
  gamma <- as.matrix(gamma_draws) %>% 
    posterior::summarise_draws() %>% 
    mutate(variable = paste0("gamma[", 1:nrow(.), "]"),
           eb_per_sec = ess_bulk/(time*60),
           et_per_sec = ess_tail/(time*60))
  
  # gamma -> ST-interaction (delta)
  if(interaction){
    delta_cats <- ncol(fit[[1]]$samples$gamma)
    for(i in 1:4){
      attr(fit[[i]]$samples$gamma, "dimnames") <- list(NULL, paste0("delta[", rep(1:delta_cats), "]"))
    }
    delta_draws <- coda::mcmc.list(fit[[1]]$samples$gamma, 
                                   fit[[2]]$samples$gamma,
                                   fit[[3]]$samples$gamma,
                                   fit[[4]]$samples$gamma)
    delta <- as.matrix(delta_draws) %>% 
      posterior::summarise_draws() %>% 
      mutate(variable = paste0("delta[", 1:nrow(.), "]"),
             eb_per_sec = ess_bulk/(time*60),
             et_per_sec = ess_tail/(time*60))
  }
  
  
  # mu
  mu <- count_draws %>% 
    posterior::summarise_draws() %>% 
    mutate(variable = paste0("mu[", 1:nrow(.), "]"),
           eb_per_sec = ess_bulk/(time*60),
           et_per_sec = ess_tail/(time*60))
  
  # Add to list
  if(interaction){
    summ <- list(beta, sigma2, rho, theta, delta, gamma, mu) %>% 
      bind_rows()
    draws <- list(beta = beta_draws,
                  sigma2 = sigma2_draws,
                  rho = rho_draws,
                  theta = theta_draws,
                  delta = delta_draws,
                  gamma = gamma_draws)
  }else{
    summ <- list(beta, sigma2, rho, theta, gamma, mu) %>% 
      bind_rows()
    draws <- list(beta = beta_draws,
                  sigma2 = sigma2_draws,
                  rho = rho_draws,
                  theta = theta_draws,
                  gamma = gamma_draws)
  }
  
  # report convergence messages
  minimum_ess <- 4 * 100
  mcmc_diags <- paste0("Median Rhat: ", round(median(summ$rhat, na.rm = T), 3), " \n",
                       round(100*mean(ifelse(summ$rhat > 1.01, 1, 0), na.rm = T), 2), "% of Rhats larger than 1.01 \n",
                       "Max Rhat = ", round(max(summ$rhat, na.rm = T), 2), 
                       " (", summ$variable[which(summ$rhat == max(summ$rhat, na.rm = T))][1], ") \n",
                       round(100*mean(ifelse(summ$ess_bulk < minimum_ess, 1, 0), na.rm = T), 2), "% of ess_bulk are too small \n", 
                       "Min ess_bulk = ", round(min(summ$ess_bulk, na.rm = T), 2), " (", 
                       summ$variable[which(summ$ess_bulk == min(summ$ess_bulk, na.rm = T))][1], ") \n",
                       round(100*mean(ifelse(summ$ess_tail < minimum_ess, 1, 0), na.rm = T), 2), "% of ess_tail are too small \n",
                       "Min ess_tail = ", round(min(summ$ess_tail, na.rm = T), 2), " (", 
                       summ$variable[which(summ$ess_tail == min(summ$ess_tail, na.rm = T))][1], ") \n",
                       "Average posterior draws per minute: ", round((n.sample*4)/time, 2))
  
  # report messages
  if(verbose) message(mcmc_diags)
  
  # MCMC settings
  nchains = 4
  MCMCsettings <- data.frame(metric = c("niter", "nburnin", "thin", "nchains", "Useable draws", 
                                        "av_pd_m", 
                                        "runTime_mins", "runTime_hours", "runTime_days"),
                             value = c(n.sample, burnin, thin, nchains, ud, 
                                       round((n.sample*nchains)/time, 2),
                                       round(time, 2), round(time/60, 2), round(time/(60*24), 2)))
  
  # Create the out object
  return(list(fit = fit,
              draws = draws,
              sres_draws = sres_draws,
              fitted_draws = fitted_draws,
              rate_draws = rate_draws,
              summary = summ,
              messages = mcmc_diags,
              modelfit = modelfit,
              MCMCsettings = MCMCsettings))
  
}

## ----------------------------------------------------------------------------
#' @param draws matrix of posterior draws for the mean of poisson
#' @param prop proportion of iterations to return a PP draw for
jf$getPoisson_yrep <- function(draws, prop = 1, seed = 45){
  
  if(!is.matrix(draws)){
    stop("Please provide a matrix with draws x n_obs!")
  }
  
  if(any(draws < 0)){
    stop("Some mu values less than zero...")
  }
  
  # take a sample of rows
  set.seed(seed)
  if(prop != 1){
	draws <- draws[sample(1:nrow(draws), size = round(prop*nrow(draws))),]
  }
  
  # get needed objects
  dims <- dim(draws)
  its <- prod(dims)
  y <- rpois(its, draws) # by columns
  
  # return the matrix of PPC
  out <- matrix(y, byrow = F, ncol = dims[2], nrow = dims[1])
  return(out)
  
}

## ----------------------------------------------------------------------------
#' @param draws matrix of posterior draws for the predicted probabilities
#' @param prop proportion of iterations to return a PP draw for
jf$getBernoulli_yrep <- function(draws, prop = 1, seed = 45){
  
  if(!is.matrix(draws)){
    stop("Please provide a matrix with draws x n_obs!")
  }
  
  if(any(draws < 0) | any(draws > 1)){
    stop("Some probabilities not in range (0, 1)...")
  }
  
  # take a sample of rows
  set.seed(seed)
  if(prop != 1){
    draws <- draws[sample(1:nrow(draws), size = round(prop*nrow(draws))),]
  }
  
  # get needed objects
  dims <- dim(draws)
  its <- prod(dims)
  y <- rbinom(its, 1, draws) # by columns
  
  # return the matrix of PPC
  out <- matrix(y, byrow = F, ncol = dims[2], nrow = dims[1])
  return(out)
  
}

## ----------------------------------------------------------------------------
#' @param y vector of observed values
#' @param yrep matrix of ppd draws (iterations by n_obs)
#' @details requires bayesplot and patchwork
jf$PoissonPPC <- function(y, yrep){

if(!any(str_detect(search(), "patchwork"))){
	stop("Please load the package 'patchwork'")
}
    
  # Create functions
  sum_y <- function(x) sum(x, na.rm = T)
  prop_zero <- function(x) mean(x == 0)
  mean_y <- function(x) mean(x, na.rm = T)
  var_y <- function(x) var(x, na.rm = T)
  
  # create plots
  sum_plot <- ppc_stat(y, yrep, stat = "sum_y") + labs(title = "Sum")
  nz_plot <- ppc_stat(y, yrep, stat = "prop_zero") + labs(title = "Proportion of\nrows with\nzero counts")
  mean_plot <- ppc_stat(y, yrep, stat = "mean_y") + labs(title = "Mean")
  var_plot <- ppc_stat(y, yrep, stat = "var_y") + labs(title = "Variance")
  
  # join plots
  (mean_plot + var_plot)/(sum_plot + nz_plot)
}

## ----------------------------------------------------------------------------
#' @details requires the `extraDistr` package
#' @param mu_draws matrix of posterior draws for the mean of poisson
#' @param p_draws matrix of posterior draws for the probability of nonzero
jf$getHurdle_yrep <- function(mu_draws, p_draws, seed = 45){

	its <- nrow(p_draws)
	yrep <- matrix(0, nrow = its, ncol = ncol(p_draws))
	
	set.seed(seed)

	pb <- txtProgressBar(min = 0, max = its, style = 3)
	for(t in 1:its){
		for(i in 1:500){
		  if(!rbernoulli(1, p_draws[t,i])){
			yrep[t,i] = 0
		  }else{
		  yrep[t,i] = extraDistr::rtpois(1, mu_draws[t,i], a = 0)
			# w = rpois(1, mu_draws[t,i])
			# while(w == 0){
			 # w = rpois(1, mu_draws[t,i]) 
			# }
			# yrep[t,i] = w
		  }
		}
		setTxtProgressBar(pb, t)
	}
	close(pb)

	# return the matrix of PPC
	return(yrep)
  
}

## ----------------------------------------------------------------------------
#' @param mu_draws matrix of posterior draws for the mean of poisson
#' @param p_draws matrix of posterior draws for the probability of nonzero
#' @param phi_draws vector of posterior draws for the overdispersion parameter
jf$getHurdleNB_yrep <- function(mu_draws, p_draws, phi_draws, seed = 45){

	its <- nrow(p_draws)
	yrep <- matrix(0, nrow = its, ncol = ncol(p_draws))
	
	set.seed(seed)
	
	pb <- txtProgressBar(min = 0, max = its, style = 3)
	for(t in 1:its){
	  for(i in 1:500){
		if(!rbernoulli(1, p_draws[t,i])){
		  yrep[t,i] = 0
		}else{
		  w = rnegbin(1, mu_draws[t,i], phi_draws[t])
		  while(w == 0){
			w = rnegbin(1, mu_draws[t,i], phi_draws[t]) 
		  }
		  yrep[t,i] = w
		}
	  }
	  setTxtProgressBar(pb, t)
	}
	close(pb)

	# return the matrix of PPC
	return(yrep)

}

## ----------------------------------------------------------------------------
#' @description provides the correct, subsetted matrix of draws from nimble run
#' @param draws all draws from nimble run (`as.matrix(fit$samples)`)
#' @param regex character expression for the param subset. For vector parameters
# use `regex = "B_qr\\["`. 
jf$getSubsetDrawsNimble <- function(draws, regex){
	draws[,which(str_detect(attr(draws, "dimnames")[[2]], regex))]
}

## ----------------------------------------------------------------------------
#' @description provides the correct, subsetted summary from 
# `SampleNimble`, `SampleCBST` and `mcmcsaeGetSummaries` 
#' @param summary summary list object from nimble or CB run (`fit$summary`)
#' @param regex character expression for the param subset. For vector parameters
# use `regex = "B_qr\\["`.
# To subset multiple parameters use `regex = "alpha|beta"` 
jf$getSubsetSummary <- function(summary, regex){
	summary %>% 
		filter(str_detect(variable, regex))
}

## ----------------------------------------------------------------------------
jf$dchi <- function(x, df, ncp = 0, log = FALSE) {
  log_f <- dchisq(x^2, df, ncp, log = TRUE) + log(2) + log(x)
  if(log) return(log_f)
  exp(log_f)
}

jf$dinvchi <- function(x, df, ncp = 0, log = FALSE) {
  log_f <- jf$dchi(1/x, df, ncp, log = TRUE) - 2*log(x)
  if(log) return(log_f)
  exp(log_f)
}

## ----------------------------------------------------------------------------
#' @param fit object from mcmcsae::MCMCMsim 
#' @param time numeric time (mins)
jf$mcmcsaeGetSummaries <- function(fit, time){
  
# subset list to just parameters of interest
  pars_draws <- fit[str_detect(names(fit), "beta|gen|linpred_")]

# function to get standard output
  foo <- function(x){
    coda::as.mcmc.list(to_mcmc(x)) %>% 
      posterior::summarise_draws()
  }

# create posterior summary object
  summ <- bind_rows(lapply(pars_draws, foo), .id = "parameter_group")
  
  # report convergence messages
  n.chains <- length(fit$beta)
  minimum_ess <- n.chains * 100
  n.sample <- dim(fit$beta[[1]])[[1]]
  mcmc_diags <- paste0("Median Rhat: ", round(median(summ$rhat, na.rm = T), 3), " \n",
                       round(100*mean(ifelse(summ$rhat > 1.01, 1, 0), na.rm = T), 2), "% of Rhats larger than 1.01 \n",
                       "Max Rhat = ", round(max(summ$rhat, na.rm = T), 2), " (", summ$variable[which(summ$rhat == max(summ$rhat, na.rm = T))][1], ") \n",
                       round(100*mean(ifelse(summ$ess_bulk < minimum_ess, 1, 0), na.rm = T), 2), "% of ess_bulk are too small \n", 
                       "Min ess_bulk = ", round(min(summ$ess_bulk, na.rm = T), 2), " (", summ$variable[which(summ$ess_bulk == min(summ$ess_bulk, na.rm = T))][1], ") \n",
                       round(100*mean(ifelse(summ$ess_tail < minimum_ess, 1, 0), na.rm = T), 2), "% of ess_tail are too small \n",
                       "Min ess_tail = ", round(min(summ$ess_tail, na.rm = T), 2), " (", summ$variable[which(summ$ess_tail == min(summ$ess_tail, na.rm = T))][1], ") \n",
                       "Average posterior draws per minute: ", round((n.sample*n.chains)/time, 2))
  # report message
  message(mcmc_diags)
  
# Output the correct object
  out <- list(summary = summ,
              messages = mcmc_diags)
  
# return object
  return(out)
  
}

## -----------------------------------------------------------------------------
#' @param summary summary list object from nimble or CB run (`fit$summary`)
#' @param regex character expression for the param subset. For vector parameters
# use `regex = "B_qr\\["`.
# To subset multiple parameters use `regex = "alpha|beta"`
#' @param n.chains integer (defaults to 4)
#' @param ess_per_chain integer (defaults to 100) required ESS per chain
jf$getSubsetConvergenceMessage <- function(summary, regex, n.chains = 4, ess_per_chain = 100){
  
  summ <- jf$getSubsetSummary(summary, regex)
  minimum_ess <- n.chains * ess_per_chain
  out <- list(medRhat = median(summ$rhat, na.rm = T),
              percRhatgr1.01 = 100*mean(ifelse(summ$rhat > 1.01, 1, 0), na.rm = T),
              percessbulkTooSmall = 100*mean(ifelse(summ$ess_bulk < minimum_ess, 1, 0), na.rm = T),
              percesstailTooSmall = 100*mean(ifelse(summ$ess_tail < minimum_ess, 1, 0), na.rm = T),
              minessbulk = min(summ$ess_bulk, na.rm = T),
              minesstail = min(summ$ess_tail, na.rm = T))
  
  # create message report
  mcmc_diags <- paste0("Number of selected parameters: ", nrow(summ), 
                       "\nNumber of chains: ", n.chains, 
                       "\nMedian Rhat: ", round(median(summ$rhat, na.rm = T), 3), " \n",
                       round(100*mean(ifelse(summ$rhat > 1.01, 1, 0), na.rm = T), 2), "% (n = ", 
                       sum(ifelse(summ$rhat > 1.01, 1, 0), na.rm = T), ") of Rhats larger than 1.01 \n",
                       "Max Rhat = ", round(max(summ$rhat, na.rm = T), 2), 
                       " (", summ$variable[which(summ$rhat == max(summ$rhat, na.rm = T))][1], ") \n",
                       round(100*mean(ifelse(summ$ess_bulk < minimum_ess, 1, 0), na.rm = T), 2), "% (n = ", 
                       sum(ifelse(summ$ess_bulk < minimum_ess, 1, 0), na.rm = T), ") of ess_bulk are too small \n", 
                       "Min ess_bulk = ", round(min(summ$ess_bulk, na.rm = T), 2), " (", 
                       summ$variable[which(summ$ess_bulk == min(summ$ess_bulk, na.rm = T))][1], ") \n",
                       round(100*mean(ifelse(summ$ess_tail < minimum_ess, 1, 0), na.rm = T), 2), "% (n = ", 
                       sum(ifelse(summ$ess_tail < minimum_ess, 1, 0), na.rm = T), ") of ess_tail are too small \n",
                       "Min ess_tail = ", round(min(summ$ess_tail, na.rm = T), 2), " (", 
                       summ$variable[which(summ$ess_tail == min(summ$ess_tail, na.rm = T))][1], ")")
  message(mcmc_diags)
  
  # output list
  return(out)
}


## ----------------------------------------------------------------------------
## Define Nimble functions
## ----------------------------------------------------------------------------

# Vectorized measurement error Poisson distribution
	dpois_me_v <- nimbleFunction(
	  run = function(x = double(1), # takes in a vector
					 mu = double(0),
					 U = double(0),
					 log = integer(0, default = 0)) {
		returnType(double(0))
		logProb <- (1/U)*sum(dpois(x, mu, log = TRUE)) # sum across the logLik
		if(log) return(logProb)
		else return(exp(logProb)) 
	  }) 
	rpois_me_v <- nimbleFunction(
	  run = function(n = integer(0, default = 1),
					 mu = double(0),
					 U = double(0)) {
		returnType(double(1))
		return(rpois(n, mu))
	  })

# vectorized poisson distribution
    dpois_v <- nimbleFunction(
		run = function(x = double(1), 
					 mu = double(1),
					 # indicator for nonzeros
					 z = double(1),
					 # Defaults to poisson
					 pos = integer(0, default = 0),
					 log = integer(0, default = 0)) {
		returnType(double(0))
		if(pos){
			# truncated poisson
			logProb <- sum( z * (dpois(x, mu, log = TRUE) - log1p(-exp(-mu))) )
		}else{
			# this is the default
			logProb <- sum( dpois(x, mu, log = TRUE) )
		}
		if(log) return(logProb)
		else return(exp(logProb)) 
	}) 
    rpois_v <- nimbleFunction(
		run = function(n = integer(0), 
					 mu = double(1),
					 z = double(1),
					 pos = integer(0, default = 0)) {
		returnType(double(1))
		return(rpois(n, mu))
	})
    
# Vectorized bernoulli
    dbern_v <- nimbleFunction(
		run = function(x = double(1), 
					 p = double(1), 
					 log = integer(0, default = 0)) {
		returnType(double(0))
		logProb <- sum(x * log(p) + (1-x) * log(1-p))
		if(log) return(logProb)
		else return(exp(logProb)) 
	}) 
    rbern_v <- nimbleFunction(
		run = function(n = integer(0), 
					 p = double(1)) {
		returnType(double(1))
		return(rbinom(n, size = rep(1, length(p)), prob = p))
	})
	  
# Vectorized pseudo-likelihood bernoulli
	dwbern_v <- nimbleFunction(
		run = function(x = double(1), 
					   w = double(1), 
					   p = double(1), 
					   log = integer(0, default = 0)) {
		  returnType(double(0))
		  logProb <- sum(w * (x * log(p) + (1-x) * log(1-p)))
		  if(log) return(logProb)
		  else return(exp(logProb)) 
	})
	rwbern_v <- nimbleFunction(
		run = function(n = integer(0), 
					   w = double(1), 
					   p = double(1)) {
		  returnType(double(1))
		  return(w * rbinom(n, size = rep(1, length(p)), prob = p))
	})
	  
# # vectorized hurdle distribution
    # dhurdle_v <- nimbleFunction(
		# run = function(x = double(1), 
					 # mu = double(1),
					 # z = double(1),
					 # p = double(1),
					 # # defaults to the hurdle distribution
					 # hurdle = integer(0, default = 1),
					 # log = integer(0, default = 0)) {
		# returnType(double(0))
		# if(hurdle){
			# # the default
			# logProb <- sum( (1-z)*log1p(-p) + z * (log(p) + dpois(x, mu, log = TRUE) - log1p(-exp(-mu))) )
		# }else{
		  # logProb <- sum( dpois(x, mu, log = TRUE) )
		# }
		# if(log) return(logProb)
		# else return(exp(logProb)) 
	# })
	# rhurdle_v <- nimbleFunction(
		# run = function(n = integer(0),
					   # mu = double(1),
					   # z = double(1),
					   # p = double(1),
					   # hurdle = integer(0, default = 1)) {
		  # returnType(double(1))
		  # return(rpois(n, mu))
	# })
		
# # vectorized NegBinom distribution
    # dnegbinom_v <- nimbleFunction(
		# run = function(x = double(1), 
					 # mu = double(1),
					 # # indicator for nonzeros
					 # z = double(1),
					 # phi = double(0),
					 # # Defaults to normal poisson
					 # pos = integer(0, default = 0),
					 # log = integer(0, default = 0)) {
		# returnType(double(0))
		# NB_lp <- lgamma(x+phi) - lfactorial(x) - lgamma(phi) + (x * (log(mu) - log(mu+phi))) + (phi * (log(phi) - log(mu+phi)))
		# NB0_lp <- phi * (log(phi) - log(mu+phi))
		# if(pos){
		  # logProb <- sum( z * ( NB_lp - log1p(-exp(NB0_lp)) ) )
		# }else{
			# # this is the default
			# logProb <- sum( NB_lp )
		# }
		# if(log) return(logProb)
		# else return(exp(logProb)) 
	# }) 
	# rnegbinom_v <- nimbleFunction(
		# run = function(n = integer(0), 
					 # mu = double(1),
					 # z = double(1),
					 # phi = double(0),
					 # pos = integer(0, default = 0)) {
		# returnType(double(1))
		# return(rpois(n, mu))
	# })