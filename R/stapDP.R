
#' Create a stapDP object
#'
#' @param object A list provided by the nd_nhpp function
#' @param chains number of chains
#' @param J number of groups
#' @return A stapDP object
#'
stapDP <- function(object,chains, J){


	out <- list(beta = coda::as.mcmc(object$beta),
				sigma = coda::as.mcmc(object$sigma))

    structure(out, class = c("stapDP"))
}
