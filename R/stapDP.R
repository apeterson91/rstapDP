
#' Create a stapDP object
#'
#' @param object A list provided by the fdp_staplm(er) functions
#' @return A stapDP object
#'
stapDP <- function(object){

	 Samples <- Parameter <- Lower <- Upper <- medp <- iteration_ix <- 
		. <- Distance <- Median <- P <-  id <- NULL
	

	K <- object$K
	spec <- object$spec
	pardf <- rbind(dplyr::tibble(iteration_ix = 1:length(object$pars$alpha),
								 Parameter = "alpha",
								 Samples = object$pars$alpha), 
			dplyr::tibble(iteration_ix = 1:length(object$pars$sigma),
						  Parameter = "sigma",
						  Samples = object$pars$sigma)) 
	probs <- object$pars$pi

	if(has_bw(spec)){
		num_penalties <- length(spec$S[[1]])
		tau_b <- object$pars$tau_b
		tau_w <- object$pars$tau_w
		colnames(tau_b) <- paste0("tau_b_",1:(K*num_penalties))
		colnames(tau_w) <- paste0("tau_w_",1:(K*num_penalties))
		scales <- cbind(tau_b,tau_w)

	}else{
		num_penalties <- length(object$spec$S)
		scales <- object$pars$scales
		colnames(scales) <- paste0("tau_",1:(K*num_penalties))
	}

	colnames(probs) <- paste0("pi","_",1:K)

	nms <- Reduce(c,lapply(spec$X,colnames))
	clnms <- Reduce(c,
					lapply(1:K,function(x) paste0("K: " , x," ",nms ))
					)

	beta <- object$pars$beta
	colnames(beta) <- c(colnames(object$mf$X),
						clnms)
	delta <- beta[,colnames(object$mf$X)]

	beta <- beta[,clnms]
	clnm_k <- lapply(1:K,function(x) grep(x=clnms,pattern = paste0("K: ",x),value=T))

	beta <- lapply(clnm_k,function(k) beta[,k])
	beta <- abind::abind(beta,along=3)



	ys <- object$pars$yhat
	gd <- expand.grid(id =paste("V_",1:ncol(object$pars$yhat)),
					  iteration_ix = 1:length(object$pars$alpha))

	colnames(ys) <- paste("V_",1:ncol(object$pars$yhat))

	ys <- suppressWarnings(dplyr::as_tibble(ys,quiet=T)) %>% 
	  dplyr::mutate(iteration_ix = 1:dplyr::n()) %>% 
		tidyr::gather(dplyr::contains("V_"),key="id",value="Samples") %>% 
		dplyr::mutate(id = as.integer(stringr::str_replace(id,"V_","")))

	yhat <- dplyr::tibble(iteration_ix = as.integer(gd$iteration_ix), 
						  Parameter = "yrep",
						  id = as.integer(gd$id))

	yhat <- suppressMessages(yhat %>% dplyr::right_join(ys))
	yhat <- rbind(yhat,dplyr::tibble(iteration_ix = rep(0,length(object$mf$y)),Parameter=rep("yobs",length(object$mf$y)), id = 1:length(object$mf$y),Samples = object$mf$y))


	out <- list(pardf = pardf,
				delta = delta,
				beta = beta,
				probs = probs,
				yhat = yhat,
				scales = scales,
				pmat = object$pars$pmat,
				cmat = object$pars$clabels,
				model = list(formula = object$formula,
							 K=(object$K),
							 y=object$mf$y,
							 alpha_a = object$alpha_a,
							 alpha_b = object$alpha_b),
				spec = object$spec,
				glmod = object$mf$glmod
				)
	if(!is.null(object$pars$subj_b)){
	  glmod <- object$mf$glmod
	  b <- object$pars$subj_b
	  D <- object$pars$subj_D
		out$subj_b <- reformat_b(b,glmod)
		out$subj_D <- reformat_D(D,glmod)
	}

    structure(out, class = c("stapDP"))
}

## ----------- internal reformatting functions
reformat_D <- function(subj_D,glmod){

	grp <- names(glmod$reTrms$cnms)
	trms <- glmod$reTrms$cnms[[1]]
	out <- subj_D

	create_nm <-function(grp,trms){
		paste0("Sigma[", grp, ": ",trms ,"]")
	}

	if(ncol(subj_D)==1){
		colnames(out) <- create_nm(grp,trms)
	}
	else if(ncol(subj_D) == 4){
		trms <- c(
		  paste0(trms[1],", ",trms[1]),
		  paste0(trms[1],", ",trms[2]),
		  paste0(trms[2],", ",trms[1]),
		  paste0(trms[2],", ",trms[2])
		  )
		colnames(out) <- create_nm(grp,trms)
	}
	return(out)
}


reformat_b <- function(subj_b,glmod){

	grp <- names(glmod$reTrms$cnms)
	trms <- sapply(glmod$reTrms$cnms[[1]],function(x) paste0(x,",", unique(glmod$reTrms$flist[[1]])))
	trms <- as.vector(trms)
	colnames(subj_b) <- paste0("b[",grp,":", trms ,"]")
	return(subj_b)
}
