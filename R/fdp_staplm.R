
#' Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @param y vector of outcomes
#' @param X design matrix
#' @param iter_ix maximum number of iterations
#' @param burn_in number of burn in iterations
#' @param thin number by which to thin samples
#' @param seed rng initializer
#' 
#' @export
#' 
fdp_staplm <- function(y,X,iter_max,burn_in,thin,seed){


	num_posterior_samples <- length(seq(from = burn_in+1, to = iter_max, by = thin))
	fit <- stapDP_backfit(y,rep(0,ncol(X)),X,iter_max,burn_in,thin,seed,num_posterior_samples)

}
