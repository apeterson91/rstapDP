
#' Retrieve Parameter Samples in Matrix Form
#'
#' Returns the parameter samples from the model associated with the 
#' fixed effects, concentration parameter, residual variance, DP probabilities
#' and smooth scales as  a S X P matrix where S is the number of samples and P is the total number of parameters
#' 
#' @param x a stapDP object
#' @param ... ignored
#' @export
#' @seealso as.array.stapDP for an array like container for the DP smooth components
#'
as.matrix.stapDP <- function(x,...){

	betamat <- Reduce(cbind,lapply(1:x$model$K,function(y){ 
									   m <- x$beta[,,y] 
									   colnames(m) <- paste0("K: ", y," ", colnames(m))
									   return(m)
						} ))
	return(cbind(x$delta,x$alpha,x$sigma,x$probs,x$scales,betamat))
}

#' Retrieve DP parameter samples in Array Form
#'
#'@export
as.array.stapDP <- function(x,...){
	return(x$beta)
}
