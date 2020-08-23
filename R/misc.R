create_S <- function(K,jg,bw,ncol_Z){

	S <- lapply(jg$pregam$S,function(x) kronecker(diag(K),x))
	S <- lapply(S,function(m) rbind(matrix(0,
										   ncol = ncol(m) + ncol_Z,
										   nrow = ncol_Z),
	           cbind(matrix(0,
							nrow=nrow(m),
							ncol=ncol_Z),
					 m))
				)
	if(bw)
		return(list(S,S))

	return(S)
}


create_X <- function(stap_term,stap_component,calc_bw,raw_X,benvo,lbls = NULL){


	X <- rbenvo::aggrenvo(benvo,raw_X,stap_term,stap_component)
	if(!is.null(lbls)){
		if(stap_component == "Distance-Time")
		  stap_component_ <- "Distance,Time"
		else
		  stap_component_ <- stap_component
		lbls_ <- stringr::str_replace(lbls,stap_component_,stap_term)
		if(calc_bw)
			lbls_ <- lapply(c("_bw","_wi"), function(x) stringr::str_replace(lbls_,stap_term,paste0(stap_term,x)))
		else{
			colnames(X) <- lbls_
			return(X)
		}
	}
	if(calc_bw){
		X <- rbenvo::bwinvo(benvo,X)
		if(!is.null(lbls))
			X <- purrr::map2(lbls_,X,function(x,y){
					colnames(y) <- x
					return(y) 
			   })
	}

	return(X)
}

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
