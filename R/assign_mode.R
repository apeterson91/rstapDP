#' Estimates Posterior Mode Cluster Assignment
#'
#' @export
#' @param x stapDP object
#' @param loss one of c("green","square") denoting which loss function should be used. Defaults to "green".
#' @param ... optional arguments passed on to specific loss function
#' @return vector of cluster assignments
#' @importFrom bendr assign_mode green_loss_unknown square_error
#'
assign_mode.stapDP <- function(x,loss = "green",...){

	if(loss=="green"){
		error <- green_loss_unknown(x$cmat,x$pmat,0.5)
		ix <- which.max(error)
	}else{
		error <- square_error(x$cmat,x$pmat)
		ix <- which.min(error)
	}
	out <- list(loss = error,
				best_loss_ix = ix,
				mode = x$cmat[ix,])

}

#' @export 
bendr::assign_mode


#' Estimates Posterior Mode Loss
#'
#' @param object stapDP object
#' @param truth true adjacency matrix
#' @param tau loss function tuning parameter
#' @param a misclassification penalty parameter
#' @param b misclassification penalty parameter - see Green and Lau for details
#' @return list of three objects (1)loss:  the loss corresponding to each clustering 
#' configuration, (2) best_loss_ix:  The index corresponding to the "best" partition,
#' and (3) mode: the best partition itself
#' @importFrom bendr green_loss
#' @export
green_loss.stapDP <- function(object,truth = NULL,tau = .5, a = 1, b = 1)
	UseMethod("green_loss")


