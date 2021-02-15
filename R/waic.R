

#' Widely Applicable Information Criterion
#' 
#' Calculates the Widely Applicable Information Criterion
#'
#' @param x model object
#' 
#' @export
waic <- function(x)
	UseMethod("waic")

#' @describeIn waic
#' @export
#' @importFrom stats dnorm
waic.stapDP <- function(x){

	warning("This function is currently experimental")
	yhatmat <- .get_yhatmat(x)
	##TODO add prior components to ll calculation
	ll <- Reduce(cbind,lapply(1:ncol(yhatmat),function(z) dnorm(x = x$model$y,
	                                                            mean = yhatmat[,z],
	                                                            sd = x$sigma[z],
	                                                            log=TRUE)))


	out <- LaplacesDemon::WAIC(ll)

	return(out)

}

.get_yhatmat <- function(x){

	iteration_ix <- Parameter <- Samples <- id <- NULL
	x$yhat %>% 
		dplyr::filter(iteration_ix>0) %>% 
		dplyr::select(-Parameter) %>% 
		tidyr::spread(iteration_ix,Samples) %>% 
		dplyr::select(-id) %>% 
		as.matrix() -> yhatmat
	
	return(yhatmat)
}
