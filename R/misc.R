
get_W <- function(glmod){

	cnms <- glmod$reTrms$cnms

	if(length(cnms)>1)
		stop("Only one grouping level currently allowed")

	cnms <- cnms[[1]]
	if(cnms[1]=="(Intercept)"){
		if(length(cnms)>1)
			return(cbind(1,glmod$fr[,cnms[2:length(cnms),drop=F]]))
		else
			return(matrix(1,nrow=nrow(glmod$fr),ncol=1))
	}
	else
		return(glmod$fr[,cnms,drop=F])

}

get_subjmat <- function(glmod){

	Matrix::fac2sparse(glmod$reTrms$flist[[1]])
}

is.mer <- function(x){
	if(!is.null(x$subj_b))
		return(TRUE)
	else
		return(FALSE)
}
