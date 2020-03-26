
#' Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @param y vector of outcomes
#' @param Z design matrix
#' @param X stap design matrix
#' @param tau_0 normal base measure scale
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param K truncation number
#' @param iter_ix maximum number of iterations
#' @param burn_in number of burn in iterations
#' @param thin number by which to thin samples
#' @param seed rng initializer
#' 
#' @export
#' 
fdp_staplm <- function(y,Z,X,
					   tau_0,
					   alpha_a,
					   alpha_b, K,
					   iter_max,burn_in,
					   thin,seed){

	stopifnot(length(tau_0) == ncol(X))

	num_posterior_samples <- sum((seq(from = burn_in+1, to = iter_max, by = 1) %% thin)==0 )
	X_ranges <- lapply(1:ncol(X),function(x) range(X[,x,drop=TRUE]))


	fit <- stapDP_fit(y,Z,X,
						  tau_0,alpha_a,alpha_b,
						  K,iter_max,burn_in,thin,
						  seed,num_posterior_samples)


	return(stapDP(fit,X_ranges))
}


#' Penalized Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @param y vector of outcomes
#' @param Z design matrix
#' @param X stap design matrix
#' @param tau_0 normal base measure scale
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param K truncation number
#' @param iter_ix maximum number of iterations
#' @param burn_in number of burn in iterations
#' @param thin number by which to thin samples
#' @param seed rng initializer
#' 
#' @export
#' 
fdp_stapplm <- function(y,Z,X,
					   alpha_a,
					   alpha_b, 
					   K,
					   nu_0,
					   iter_max,burn_in,
					   thin,seed){

	stopifnot(nu_0>0)

  num_posterior_samples <- sum((seq(from = burn_in+1, to = iter_max, by = 1) %% thin)==0 )
  X_ranges <- lapply(1:ncol(X),function(x) range(X[,x,drop=TRUE]))


	fit <- stappDP_fit(y,Z,X,
					   nu_0,alpha_a,alpha_b,
					   K,iter_max,burn_in,thin,
					   seed,num_posterior_samples)


	return(stapDP(fit,X_ranges))
}
