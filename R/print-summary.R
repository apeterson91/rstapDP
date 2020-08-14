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
	cat("fdp_staplm \n")
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
	cat("\n Auxiliary Variables \n ")
	
	mat <- x$pardf %>% tidyr::spread(Parameter,Samples) %>% dplyr::select(-iteration_ix) %>% 
		as.matrix()

	.printfr(.median_and_madsd(mat),digits)

	if(is.mer(x)){
      cat("\nError terms:\n")
	  foo <- VarCorr(x)
	  print(foo)
	  cat("Num. levels:",
          paste(names(ngrps(x)), unname(ngrps(x)), collapse = ", "), "\n")
	}


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

