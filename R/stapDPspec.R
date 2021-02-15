#  This software is part of the rstapDP package
#  Copyright (C) 2020 Adam Peterson
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Retrieve STAP-DP model specification from STAP model formula
#'
#' Get stapDPspec object which details the various components of the stapDP model specification
#'
#' @param f formula from \code{\link{fdp_staplm}}, \code{\link{fdp_staplmer}}
#' @param K DP truncation integer
#' @param benvo Built Environment object - \code{\link[rbenvo]{benvo}} - containing data for model 
#' @return \code{\link{stapDPspec}} object
#'
get_stapDPspec <- function(f,K,benvo){

    with_bars <- lme4::findbars(f)
    f <- lme4::nobars(f)
	get_ics <- function(f,vec_var){
		which(all.names(f) %in% vec_var)
	}
	get_k <- function(strings){
	  out <- stringr::str_extract(strings,", ?k ?= ?[1-9][1-9]? *\\)")
	  out <- sapply(out,function(x) if(is.na(x)) return(")") else x)
	  return(out)
	}
	get_names <- function(f,ics){
		all.names(f)[ics +1] 
	}
	get_indicator <- function(f,ics,vec_var){
		(all.names(f)[ics] %in% vec_var)*1
	}
    stap_ics <- get_ics(f, c("stap","stap_bw"))
    sap_ics <- get_ics(f,c("sap","sap_bw"))
    tap_ics <- get_ics(f,c("tap","tap_bw"))
    if(!length(stap_ics) & !length(sap_ics) & !length(tap_ics))
        stop("No covariates designated as 'stap','sap',or 'tap'  in formula", .call = F)
	stap_nms <- get_names(f,stap_ics)
	stap_bw <- get_indicator(f,stap_ics, c("stap_bw"))
	sap_nms <- get_names(f,sap_ics)
	sap_bw <- get_indicator(f,sap_ics,c("sap_bw"))
    tap_nms <- get_names(f,tap_ics)
	tap_bw <- get_indicator(f,tap_ics,c("tap_bw")) 
	tms <- attr(terms(f),"term.labels")

	stap_tms <- tms[stringr::str_detect(tms,"(^stap\\()|(^stap_bw\\()")]
	tap_tms <- tms[stringr::str_detect(tms,"(^tap\\()|(^tap_bw\\()")]
	sap_tms <- tms[stringr::str_detect(tms,"(^sap\\()|(^sap_bw\\()")]

	stap_k <- get_k(stap_tms)
	tap_k <- get_k(tap_tms)
	sap_k <- get_k(sap_tms)
	
	if(length(stap_nms)>0){
		stap_nms <- cbind(stap_nms,"Distance-Time",stap_bw,stap_k)
	}
	if(length(sap_nms)>0)
		sap_nms <- cbind(sap_nms,"Distance",sap_bw,sap_k)
	if(length(tap_nms)>0)
		tap_nms <- cbind(tap_nms,"Time",tap_bw,tap_k)

	stap_mat <-rbind(stap_nms,sap_nms,tap_nms)

    not_needed <- c(stap_nms,sap_nms,tap_nms)
    formula_components <- all.vars(f)[!(all.vars(f) %in% not_needed)]
    if(!attr(terms(f),"intercept"))
        formula_components <- c(formula_components,"0")
    if(grepl("cbind",all.names(f))[2]){
        new_f1 <- paste0("cbind(",formula_components[1],", ",formula_components[2], ")", " ~ ")
        ix <- 3
    }
    else{
        new_f1 <- paste0(formula_components[1],' ~ ')
        ix <- 2
    }

    new_f2 <- paste(formula_components[ix:length(formula_components)],collapse = "+")
    new_f <- paste0(new_f1,new_f2)
	if(length(with_bars)){
		mer_f <- paste0(lapply(with_bars,function(x) paste0("(",deparse(x),")")),collapse = " + ")
		new_f <- paste0(new_f," + ",mer_f)
	}

    return(
		   stapDPspec(stapless_formula = as.formula(new_f, env = environment(f)),
					  stap_mat = stap_mat,
					  K = K,
					  benvo = benvo
			   )
		   )
}


#' Create STAP-DP data structure
#' 
#' @param stapless_formula from \code{\link{get_stapDPspec}}
#' @param stap_mat matrix of stap specification properties 
#' @param K DP truncation
#' @param benvo Built Environment object - \code{\link[rbenvo]{benvo}} - containing data for model 
#'
stapDPspec <- function(stapless_formula,stap_mat,K,benvo){


	term <- stap_mat[,1]
	component <- stap_mat[,2]
	between_within <- as.integer(stap_mat[,3])
	dimension <- sapply(stap_mat[,4],function(x){ 
							if(stringr::str_length(x)==1) 
								return(10)
							else 
								as.integer(stringr::str_replace(stringr::str_replace(stap_mat[,4],
																					 "\\)",""),
																", ?k = ",""))
		   })

	if(!(all(unique(term)==term)))
		stop("Only one BEF name may be assigned to a stap term e.g. no sap(foo) + tap(foo)\n
			 If you wish to model components this way create a different name e.g. sap(foo) + tap(foo_bar)")
	if(!all(term %in% rbenvo::bef_names(benvo)))
		stop("All stap terms must have data with corresponding name in benvo")
	if(length(term)>1)
		stop("Only one stap/sap/tap term allowed")
	## adapting for iid gaussian priors

create_unique_ID_mat <- function(id_one,id_two = NULL){
	tmp <- paste0(id_one,"_",id_two)
	lvls <- unique(tmp)
	new_id <- factor(tmp,levels=lvls)
	Matrix::fac2sparse(new_id)
}

	temp_df <- rbenvo::joinvo(benvo,term = term,
							  component = component,
							  NA_to_zero = FALSE)
	ids <- rbenvo::get_id(benvo)

	Xmat <- temp_df %>% 
		dplyr::group_by_at(ids) %>% 
		dplyr::mutate(bef_id = stringr::str_c("BEF_",1:dplyr::n())) %>%
		tidyr::pivot_wider(id_cols = ids,
						   names_from = "bef_id" ,
						   values_from = component) %>% 
	  dplyr::ungroup() %>% 
		dplyr::select_at(dplyr::vars(!dplyr::contains(ids))) %>%
		as.matrix()

	if(length(ids)>1)
		idmat <- create_unique_ID_mat(benvo$subject_data[,ids[1],drop=T],
									  benvo$subject_data[,ids[2],drop=T])
	else
		idmat <- Matrix::fac2sparse(benvo$subject_data[,ids,drop=T])

	Xmat <- as.matrix(Matrix::t(idmat) %*% Xmat)
	msk <- which(is.na(Xmat))
	Xmat[msk] <- -1
	L <- (Xmat> 0)*1

	
	noise <- stats::rnorm(nrow(Xmat))
  
	out <- mgcv::jagam(formula = noise ~ 0 + s(Xmat,by = L,bs='ps',k = dimension), family = gaussian(), 
					   data = temp_df,
					   file = tempfile(fileext = ".jags"), 
					   offset = NULL,
					   centred = FALSE,
					   diagonalize = TRUE)
	ranges <- list()
	ranges[[term]] <- list()
	if(component=="Distance"|component=="Distance-Time")
		ranges[[term]]$Distance <- range(benvo$sub_bef_data[[term]]$Distance,na.rm=T)
	if(component=="Time"|component=="Distance-Time")
		ranges[[term]]$Time = range(benvo$sub_bef_data[[term]]$Time,na.rm=T)

	mf <- rbenvo::subject_design(benvo,lme4::nobars(stapless_formula))

	f_ <-  lme4::nobars(stapless_formula)
	num_fixed <- ncol(mf$X)

	X <- out$jags.data$X
	nms <- stringr::str_c("s(",term,".",1:ncol(X),")")
	colnames(X) <- nms
	if(between_within)
		X <- .create_X_bw_wi(X,
		                     Matrix::fac2sparse(benvo$subject_data[,ids[1],drop=T]),
		                     term)
	
	



	sobj <- out$pregam$smooth
	sobj[[1]]$term <- component

	out <- list(stapless_formula = stapless_formula,
				term = term,
				component = component,
				between_within = between_within,
				dimension = dimension,
				ranges = ranges,
				X = X,
				mf = mf,
				smooth_objs = sobj
				)

	structure(out,class=c("stapDPspec"))
}


#------------------------


.create_X_bw_wi <- function(X, idmat,term){

	Xb <- Matrix::t(idmat) %*% (idmat %*% X)
	Xw <- X - Xb
	Xb <- as.matrix(Xb)
	Xw <- as.matrix(Xw)
	colnames(Xb) <- stringr::str_c("s(",term,"_bw.",1:ncol(Xb),")")
	colnames(Xw) <- stringr::str_c("s(",term,"_wi.",1:ncol(Xw),")")


	return(list(X_bw = Xb,X_wi=Xw))


}
