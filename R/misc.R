
#' STAP Smooth Term
#' 
#' @export
#' @param term single covariate name to be clustered 
#' @param k basis dimension see \code{\link[mgcv]{s}}.
#'
s_stap <- function(term,  k){
	foo <- mgcv::s(term, k = k,bs='ps',fx=FALSE,
				   m=c(2,2),by=NA,xt=NULL,
				   id=NULL,sp=NULL,pc=NULL)

	foo$term <- term
	foo$label <- paste0("s(", term, ")")
	return(foo)
}
