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
	cat("\n Observations:", nrow(x$pmat))
	cat("\n Formula: ", formula_string(x$model$formula))
	cat("\n Fixed predictors: ",length(unique(x$fixef$Parameter)))
	cat("\n K: ",as.character(x$model$K))
	
	cat("\n ----------------------- \n ")

	cat("\n Fixed Effects  \n ")

	mat<- x$fixef %>% tidyr::spread(Parameter,Samples)  %>% dplyr::select(-iteration_ix) %>% 
		as.matrix()

	.printfr(t(.median_and_madsd(mat)),digits)

	cat("\n ----------------------- \n ")

	cat("Number of Clusters\n ")

	cat(summary(apply(x$cmat,1,function(x) length(unique(x)))))

	topKs <- x$probs %>% dplyr::group_by(K) %>% dplyr::summarise(med = median(Samples)) %>% dplyr::filter(med>0.001) %>% dplyr::arrange(dplyr::desc(med)) %>% utils::head(5) %>% 
		dplyr::pull(K)

	cat("\n ----------------------- \n ")

	cat("\n Cluster Probabilities  \n ")

	.printfr(t(.median_and_madsd(x$probs %>% dplyr::select(-Parameter) %>% dplyr::filter(K %in% topKs) %>%
								 tidyr::spread(K,Samples) %>% dplyr::select(-iteration_ix))),digits)


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
#' @keywords internal
#' @param x a stapDP object
#'
diagnostics <- function(x)
	UseMethod("diagnostics")

#'
#' 
#' @export
#' @describeIn diagnostics
#' 
diagnostics.stapDP <- function(x){

	Parameter <- Samples<- iteration_ix <- select <- NULL
	alphasamps <- x$pardf %>% dplyr::filter(Parameter == 'alpha' ) %>% dplyr::pull() 
	alpha <- coda::as.mcmc(alphasamps)

	sigma <- coda::as.mcmc(x$pardf %>% dplyr::filter(Parameter == 'sigma' ) %>% dplyr::pull())

	K <- max(as.integer(x$probs$K))

	lbls <- paste0("pi_",1:K)
	pis <- coda::as.mcmc(x$probs %>% tidyr::spread(Parameter,Samples)  %>%
						 dplyr::mutate_at(lbls,function(x) tidyr::replace_na(x,0)) %>%
						 dplyr::group_by(iteration_ix) %>% 
						 dplyr::summarise_at(lbls,sum) %>% 
						   select(lbls) %>% as.matrix())

	matone <- c(coda::geweke.diag(alpha)[1]$z,coda::geweke.diag(sigma)[1]$z)
	names(matone) <- c("alpha","sigma")
	mattwo <- coda::geweke.diag(pis)[1]$z

	cat("Diagnostics: \n ")
	if(all(alphasamps==0)){
		cat("Alpha has collapsed at zero. Expect to see NaN values for pis~=zero ")
		cat("\n")
	}
	print(matone)
	cat("\n")
	print(mattwo)
	return(invisible(list(alpha_sigma = matone,pis = mattwo)))
}
