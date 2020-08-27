#' Methods for stapDP objects
#'
#' @templateVar stapDPregArg object, x 
#' @template args-stapDP-object
#' @template args-dots-ignored 
#'
#' @title Methods for stapDP objects
#' @description Methods for stapDP objects that are similar to
#' their counterparts in the \pkg{stats} or \pkg{lme4} packages.
#' @details Almost all methods behave intuitively as one familiar with the \pkg{stats} or \pkg{lme4} 
#' packages would expect. 
#' @name stapDP-methods
#' @aliases VarCorr ngrps sigma nsamples
#'
#'
#' @importFrom stats coef nobs formula
#' 
NULL

#' @rdname stapDP-methods
#' @export
coef.stapDP <- function(object, ...) {
 apply(as.matrix(object),2,median)
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
	sc <- x$sigma
	Sigma_ <- colMeans(x$subj_D[,grepl("^Sigma\\[", colnames(x$subj_D)), drop = FALSE])
	nc <- vapply(cnms, FUN = length, FUN.VALUE = 1L)
	nms <- names(cnms)
	ncseq <- seq_along(nc)
	Sigma <- matrix(Sigma_,nc,nc)
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


#' @rdname stapDP-methods
#' @export
#' @importFrom lme4 ranef
#' @export ranef
#'
ranef.stapDP <- function(object,...){

	.glmer_check(object)
	b_nms <- colnames(object$subj_b)
	point_estimates <- summary(object)[b_nms,"50%"]
	out <- ranef_template(object)
	group_vars <- names(out)
	for (j in seq_along(out)) {
		tmp <- out[[j]]
		pars <- colnames(tmp) 
		levs <- rownames(tmp)
		levs <- gsub(" ", "_", levs) 
	for (p in seq_along(pars)) {
	  stan_pars <- paste0("b[",group_vars[j], ":", pars[p],",", levs, "]")
	  tmp[[pars[p]]] <- unname(point_estimates[stan_pars])
	}
		out[[j]] <- tmp
	}
	out
}

# Call lme4 to get the right structure for ranef objects
#' @importFrom lme4 lmerControl glmerControl lmer glmer 
ranef_template <- function(object) {
  
    new_formula <- object$spec$stapless_formula 
	lme4_fun <- "lmer"

	cntrl_args <- list(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1))
	cntrl_args$check.conv.grad <- "ignore"
	cntrl_args$check.conv.singular <- "ignore"
	cntrl_args$check.conv.hess <- "ignore"
	cntrl_args$check.nlev.gtreq.5 <- "ignore"
	cntrl_args$check.nobs.vs.rankZ <- "ignore"
	cntrl_args$check.nobs.vs.nlev <- "ignore"
	cntrl_args$check.nobs.vs.nRE <- "ignore"
  
	cntrl <- do.call(paste0(lme4_fun, "Control"), cntrl_args)
  
	fit_args <- list(
	formula = new_formula,
	data = object$model$benvo@subject_data,
	control = cntrl
	)
  
	lme4_fit <- suppressWarnings(do.call(lme4_fun, args = fit_args))
	ranef(lme4_fit)
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
