#' Estimates Posterior Mode Cluster Assignment
#'
#' @export
#' @param x stapDP object
#' @param loss one of c("green","square") denoting which loss function should
#' be used. Defaults to "green".
#' @param ... optional arguments passed on to specific loss function
#' @return vector of cluster assignments
#' @importFrom bendr assign_mode green_loss_unknown square_error
#'
assign_mode.stapDP <- function(x, loss = "green", ...) {
  if (loss == "green") {
    error <- green_loss_unknown(x$cmat, x$pmat, 0.5)
    ix <- which.max(error)
  } else {
    error <- square_error(x$cmat, x$pmat)
    ix <- which.min(error)
  }
  out <- list(loss = error, best_loss_ix = ix, mode = x$cmat[ix, ])
  return(out)
}

#' @export
bendr::assign_mode


#' @describeIn green_loss
#' @importFrom bendr green_loss_known
#' @export
green_loss.stapDP <- function(object, truth = NULL, tau = 0.5, a = 1, b = 1) {
  if (is.null(truth)) {
    loss <- green_loss_unknown(object$cmat, object$pmat, tau)
    ix <- which.max(loss)
  } else {
    stopifnot(dim(truth)[1] == ncol(object$cmat))
    stopifnot(dim(truth)[2] == ncol(object$cmat))
    loss <- green_loss_known(object$cmat, object$pmat, truth, a, b)
    ix <- which.min(loss)
  }
  out <- list(loss = loss, best_loss_ix = ix, mode = object$cmat[ix, ])
  return(out)
}


#' @export
#' @importFrom bendr green_loss
bendr::green_loss
