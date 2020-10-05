#' Print and Summary Method
#'
#' @export
#' @method print stapDP
#' @importFrom stats quantile mad median quantile
#' @param x stapDP object
#' @param digits number of digits to round off to
#' @param ...  ignored
#' 
print.stapDP <- function(x,digits=1,...){


	Parameter <- Samples <- iteration_ix <- K <- P <- med <- NULL
	if(is.mer(x))
		cat("fdp_staplmer")
	else
		cat("fdp_staplm")
	cat("\n Observations:", length(x$model$y))
	cat("\n Formula: ", formula_string(x$model$formula))
	cat("\n Fixed predictors: ", ncol(x$delta))
	cat("\n K: ",as.character(x$model$K))
	
	cat("\n ----------------------- \n ")

	cat("\n Fixed Effects  \n ")

	mat<- x$delta

	.printfr(t(.median_and_madsd(mat)),digits)

	cat("\n ----------------------- \n ")

	cat("Number of Clusters\n ")
	cnumsum <- summary(apply(x$cmat,1,function(x) length(unique(x))))
	cat(names(cnumsum))
	cat("\n")
	cat(cnumsum)

	cat("\n ----------------------- \n ")

	cat("\n Cluster Probabilities  \n ")

	.printfr(t(.median_and_madsd(x$probs)) ,digits)


	cat("\n ----------------------- \n ")
	cat("\n Concentration Parameter \n ")
	

	.printfr(.median_and_madsd(x$alpha),digits)

	if(is.mer(x)){
      cat("\nError terms:\n")
	  foo <- VarCorr(x)
	  print(foo)
	  cat("Num. levels:",
          paste(names(ngrps(x)), unname(ngrps(x)), collapse = ", "), "\n")
	}else{
      cat("\nResidual sd:\n")
		.printfr(.median_and_madsd(x$sigma),digits)
	}


}


#' Summary method for stapDP objects
#' 
#' Summaries of parameter estimates and MCMC convergence diagnostics 
#' (effective sample size, Rhat).
#'
#' @export
#' @method summary stapDP
#'
#' @param object stapDP object
#' @param ... ignored
#' @param digits Number of digits to use for formatting numbers when printing. 
#'   When calling \code{summary}, the value of digits is stored as the 
#'   \code{"print.digits"} attribute of the returned object.
#' @return The \code{summary} method returns an object of class 
#'   \code{"summary.staDP"}  which is a matrix of 
#'   summary statistics and diagnostics, with attributes storing information for use by the
#'   \code{print} method. 
#'
summary.stapDP <- function(object,...,digits=1){


	out <- object$summary
    stats <- colnames(out)
    if ("n_eff" %in% stats) {
      out[, "n_eff"] <- round(out[, "n_eff"])
    }

	structure(
	out,
	call = object$call,
	family = "gaussian[identity]",
	formula = formula(object),
    posterior_sample_size = nsamples(object),
    nobs = nobs(object),
	ngrps = if(is.mer(object)) ngrps(object) else NULL,
	print.digits = digits,
	class = "summary.stapDP"
	)
}

#' @rdname summary.stapDP
#' @export
#' @method print summary.stapDP
#'
#' @param x An object of class \code{"summary.stapDP"}.
print.summary.stapDP <-
  function(x, digits = max(1, attr(x, "print.digits")),
           ...) {

    atts <- attributes(x)
    cat("\nModel Info:")
    cat("\n family:      ", atts$family)
    cat("\n formula:     ", formula_string(atts$formula))
	cat("\n sample:      ", atts$posterior_sample_size, 
	  "(posterior sample size)")
    
    cat("\n observations:", atts$nobs)
    if (!is.null(atts$ngrps)) {
      cat("\n groups:      ", paste0(names(atts$ngrps), " (", 
                                     unname(atts$ngrps), ")", 
                                     collapse = ", "))
    }
	cat("\n")
    
	hat <- "Rhat"
	str_diag <- "MCMC diagnostics"
	str1 <- "and Rhat is the potential scale reduction factor on split chains"
	str2 <- " (at convergence Rhat=1).\n"
    sel <- which(colnames(x) %in% c("mcse", "n_eff", hat))
    has_mc_diagnostic <- length(sel) > 0
    if (has_mc_diagnostic) {
      xtemp <- x[, -sel, drop = FALSE]
      colnames(xtemp) <- paste(" ", colnames(xtemp))
    } else {
      xtemp <- x
    }
    
    # print table of parameter stats
    .printfr(xtemp, digits)
    
    if (has_mc_diagnostic) {
      cat("\n", str_diag, "\n", sep = '')
      mcse_hat <- format(round(x[, c(hat), drop = FALSE], digits), 
                          nsmall = digits)
      n_eff <- format(x[, "n_eff", drop = FALSE], drop0trailing = TRUE)
      print(cbind(mcse_hat, n_eff), quote = FALSE)
      cat("\n n_eff is a crude measure of effective sample size, ", 
          str1, 
          str2, sep = '')
    }
    
    invisible(x)
  }

# internal use -------------------------------------------------------------------------------------

.printfr <- function(x, digits, ...) {
  print(format(round(x, digits), nsmall = digits), quote = FALSE, ...)
}
.median_and_madsd <- function(x) {
  cbind(Median = apply(x, 2, median), MAD_SD = apply(x, 2, mad))
}

formula_string <- function(formula, break_and_indent = TRUE) {
  coll <- if (break_and_indent) "--MARK--" else " "
  char <- gsub("\\s+", " ", paste(deparse(formula), collapse = coll))
  if (!break_and_indent)
    return(char)
  gsub("--MARK--", "\n\t  ", char, fixed = TRUE)
}

