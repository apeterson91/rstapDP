#' Methods for stapDP objects
#'
#'
#' @name stapDP-methods
#' @aliases VarCorr ngrps sigma nsamples
#'
#' @param object stapDP object
#' @param ... Ignored
#'
#' @importFrom stats coef nobs formula
#' 
NULL

#' @rdname stapDP-methods
#' @export
coef.stapDP <- function(object, ...) {
 apply(object$delta,2,median)
}

#' @rdname stapDP-methods
#' @export
#'
confint.stapDP <- function(object, ...) {
    stop("Please use posterior_interval() to obtain ", 
         "Bayesian interval estimates.", 
         call. = FALSE)
}

#' @rdname stapDP-methods
#' @export
fitted.stapDP <- function(object, ...)  
    return(object$yhat)


#' @rdname stapDP-methods
#'
#' @export
nobs.stapDP <- function(object, ...){
	length(object$model$y)
}

#' @rdname stapDP-methods
#'
#' @param x stapDP object
#' @export
#' @importFrom nlme VarCorr
#' @importFrom stats cov2cor
#'
VarCorr.stapDP <- function(x){

	cnms <- .cnms(x)
	sc <- x$pardf %>% dplyr::filter(Parameter == 'sigma') %>% dplyr::pull(Samples)
	Sigma_ <- colMeans(x$subj_D[,grepl("^Sigma\\[", colnames(x$subj_D)), drop = FALSE])
	nc <- vapply(cnms, FUN = length, FUN.VALUE = 1L)
	nms <- names(cnms)
	ncseq <- seq_along(nc)
	Sigma <- matrix(Sigma_,1,1)
	stddev <- sqrt(diag(Sigma))
	corr <- cov2cor(Sigma)
	ans <- list(structure(Sigma , stddev = stddev,correlation = corr))
	names(ans) <- nms

	structure(ans, sc = mean(sc), useSc = TRUE, class = "VarCorr.merMod")

}

#' Formula method for stapDP objects
#'
#' @keywords internal
#' @export
#' @param x a stapDP object
#' @param ... ignored currently
formula.stapDP <- function(x,...){
	x$formula
}

#' @rdname stapDP-methods
#'
#' @export
#' @export ngrps
#' @importFrom lme4 ngrps
#' 
ngrps.stapDP <- function(object, ...) {
  vapply(.flist(object), nlevels, 1)  
}


#' @rdname stapDP-methods
#' @export
#' @importFrom rstantools nsamples
nsamples.stapDP <- function(object, ...) {
  nrow(object$delta)
}





.cnms <- function(object, ...) UseMethod(".cnms")
.cnms.stapDP <- function(object, ...) {
  .glmer_check(object)
  object$glmod$reTrms$cnms
}

.flist <- function(object, ...) UseMethod(".flist")
.flist.stapDP <- function(object, ...) {
  .glmer_check(object)
  as.list(object$glmod$reTrms$flist)
}

.glmer_check <- function(object) {
  if (!is.mer(object))
    stop("This method is for fdp_staplmer models only.", 
         call. = FALSE)
}
