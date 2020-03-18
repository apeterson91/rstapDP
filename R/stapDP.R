
#' Create a stapDP object
#'
#' @param object A list provided by the nd_nhpp function
#' @param chains number of chains
#' @param J number of groups
#' @return A stapDP object
#'
stapDP <- function(object,X_ranges){


	out <- list(beta = coda::as.mcmc(object$beta),
				pi = coda::as.mcmc(object$pi),
				sigma = coda::as.mcmc(object$sigma),
				alpha = coda::as.mcmc(object$alpha),
				cluster_assignment = coda::as.mcmc(object$cluster_assignment),
				pmat = object$PairwiseProbabilityMat,
				X_ranges = X_ranges
	)

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
