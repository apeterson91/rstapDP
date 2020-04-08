
#' Create a stapDP object
#'
#' @param object A list provided by the fdp_stap(p)lm function
#' @param X_ranges the range of the covariate space for the STAP
#' @return A stapDP object
#'
stapDP <- function(object,X_ranges){


	beta <- object$beta
	colnames(beta) <- paste0("beta_",1:ncol(beta))
	pi <- object$pi
	colnames(pi) <- paste0("pi_",1:ncol(pi))
	out <- list(beta = coda::as.mcmc(beta),
				pi = coda::as.mcmc(object$pi),
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
	if(!is.null(object$sigma)){
		sigma <- object$sigma
		names(sigma) <- "sigma"
		out$sigma <- coda::as.mcmc(sigma)
	}

    structure(out, class = c("stapDP"))
}

