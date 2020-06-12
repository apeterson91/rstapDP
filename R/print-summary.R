#' Print and Summary Method
#'
#' @export
#' @importFrom stats quantile mad median quantile
#' @param x stapDP object
#' @param digits number of digits to round off to
#' 
print.stapDP <- function(x,digits=1){


	Parameter <- Samples <- iteration_ix <- K <- P <- NULL
	cat("fdp_staplm \n")
	cat("\n Observations:", nrow(x$pmat))
	cat("\n Formula: ", formula_string(x$model$formula))
	cat("\n Fixed predictors: ",length(unique(x$fixef$Parameter)))
	cat("\n K: ",as.character(x$model$K))
	
	cat("\n ----------------------- \n ")

	cat("\n Fixed Effects  \n ")

	mat<- x$fixef %>% tidyr::spread(Parameter,Samples)  %>% dplyr::select(-iteration_ix) %>% 
		as.matrix()

	.printfr(t(.median_and_madsd(mat)),digits)

	mat <- x$pardf %>% tidyr::spread(Parameter,Samples) %>% dplyr::select(-iteration_ix) %>% 
		as.matrix()

	.printfr(.median_and_madsd(mat),digits)

	cat("\n ----------------------- \n ")

	cat("\n Random Effects  \n ")

	mat<- x$ranef %>% dplyr::select(-K,-P) %>% 
		tidyr::spread(Parameter,Samples)  %>% dplyr::select(-iteration_ix) %>% 
		as.matrix()
	for(k in 0:(x$model$K-1)){
  		p <- ncol(mat) / x$model$K
  		start <- (k*p+1)
  		end <- start + p - 1
  	  .printfr(t(.median_and_madsd(mat[,start:end,drop=F])),digits)
		}

	cat("\n Scales  \n ")

	.printfr(t(.median_and_madsd(x$scales)),digits)


	cat("\n Cluster Probabilities  \n ")

	.printfr(t(.median_and_madsd(x$probs %>% dplyr::select(-K) %>% tidyr::spread(Parameter,Samples) %>% dplyr::select(-iteration_ix))),digits)

	cat("\n ----------------------- \n ")
	cat("\n Auxiliary Variables \n ")
	
	mat <- x$pardf %>% tidyr::spread(Parameter,Samples) %>% dplyr::select(-iteration_ix) %>% 
		as.matrix()

	.printfr(.median_and_madsd(mat),digits)
	

}

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

#' Diagnostics
#'
#' @export
#' 
diagnostics <- function(x)
	UseMethod("diagnostics")

#' Diagnostics 
#'
#' @describeIn diagnostics
#' 
diagnostics.stapDP <- function(x){

	Parameter <- NULL

	alpha <- coda::as.mcmc(x$pardf %>% dplyr::filter(Parameter == 'alpha' ) %>% dplyr::pull())

	sigma <- coda::as.mcmc(x$pardf %>% dplyr::filter(Parameter == 'sigma' ) %>% dplyr::pull())

	mat <- c(coda::geweke.diag(alpha)[1]$z,coda::geweke.diag(sigma)[1]$z)
	names(mat) <- c("alpha","sigma")

	cat("Diagnostics: \n ")
	print(mat)
}
