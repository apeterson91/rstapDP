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



#' @export
#' @importFrom bendr green_loss
bendr::green_loss
