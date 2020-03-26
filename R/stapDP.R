
#' Create a stapDP object
#'
#' @param object A list provided by the nd_nhpp function
#' @param chains number of chains
#' @param J number of groups
#' @return A stapDP object
#'
stapDP <- function(object,X_ranges){


	beta <- object$beta
	colnames(beta) <- paste0("beta_",1:ncol(beta))
	pi <- object$pi
	colnames(pi) <- paste0("pi_",1:ncol(pi))
	out <- list(beta = coda::as.mcmc(beta),
				pi = coda::as.mcmc(object$pi),
				sigma = coda::as.mcmc(object$sigma),
				alpha = coda::as.mcmc(object$alpha),
				cluster_assignment = coda::as.mcmc(object$cluster_assignment),
				pmat = object$PairwiseProbabilityMat,
				X_ranges = X_ranges
	)
	if(!is.null(object$tau)){
		tau <- object$tau
		colnames(tau) <- paste0("tau_",1:ncol(tau))
		out$tau <- coda::as.mcmc(tau)
	}

    structure(out, class = c("stapDP"))
}


# get_par_tibble <- function(object) UseMethod("get_par_tibble")
# 
# get_par_tibble.stapDP <- function(object){
#   
#   tibble(Iteration = 1:nrow(object$beta),
#          Parameter = )
#   
# }
