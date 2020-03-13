
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
					   tau_0,alpha_a,
					   alpha_b, K,
					   iter_max,burn_in,
					   scale_center = TRUE,
					   thin,seed){

	stopifnot(length(tau_0) == ncol(X))

	num_posterior_samples <- length(seq(from = burn_in+1, to = iter_max, by = thin))
	X_ranges <- lapply(1:ncol(X),function(x) range(X[,x,drop=TRUE]))

	if(scale_center){
  	xbar <- colMeans(X)
  	zbar <- mean(Z[,1])
  	xsd <- apply(X,2,sd)
  	scaled_Z <- cbind(Z[,1],Z[,2] - zbar)
  	scaled_X <- (X - xbar) / xsd
	}
	else{
	  scaled_Z <- Z
	  scaled_X <- X
	}
	


	fit <- stapDP_backfit(y,scaled_Z,scaled_X,
						  tau_0,alpha_a,alpha_b,
						  K,iter_max,burn_in,thin,
						  seed,num_posterior_samples)

	if(scale_center){
	  fit$beta[,1] <- fit$beta[,1] - fit$beta[,2] * zbar 
	  
  	for(i in 3:(ncol(fit$beta)-2))
  	  fit$beta[,1] <- fit$beta[,1] - rowSums(fit$beta[,i:(i+2)] * xbar) / xsd
  	
  	for(i in 3:(ncol(fit$beta)-2))
  	  fit$beta[,i:(i+2)] <- fit$beta[,i:(i+2)] / xsd 
	}
	
	

	return(stapDP(fit,X_ranges))
}
