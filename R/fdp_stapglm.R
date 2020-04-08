#' Penalized Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Logistic Linear Model
#' 
#' @param y vector of outcomes
#' @param n vector of trials
#' @param Z design matrix
#' @param X stap design matrix
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param K truncation number
#' @param nu_0 df for coefficients scale prior
#' @param iter_max maximum number of iterations
#' @param burn_in number of burn in iterations
#' @param thin number by which to thin samples
#' @param seed rng initializer
#' 
fdp_stappglm <- function(y,n,Z,X,
					   alpha_a,
					   alpha_b, 
					   K,
					   nu_0,
					   iter_max,burn_in,
					   thin,seed){

	stopifnot(nu_0>0)

  num_posterior_samples <- sum((seq(from = burn_in+1, to = iter_max, by = 1) %% thin)==0 )
  X_ranges <- lapply(1:ncol(X),function(x) range(X[,x,drop=TRUE]))


#	fit <- stappDP_logistic_fit(y,n,
#								Z,X,
#							   nu_0,alpha_a,alpha_b,
#							   K,iter_max,burn_in,thin,
#							   seed,num_posterior_samples)
  fit <- NULL


	return(stapDP(fit,X_ranges))
}
