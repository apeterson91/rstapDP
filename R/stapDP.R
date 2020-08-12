
#' Create a stapDP object
#'
#' @param object A list provided by the fdp_staplm.fit function
#' @return A stapDP object
#'
stapDP <- function(object){

	 Samples <- Parameter <- Lower <- Upper <- medp <- iteration_ix <- 
		. <- Distance <- Median <- P <-  id <- NULL
	

	K <- object$K
	pardf <- rbind(dplyr::tibble(iteration_ix = 1:length(object$pars$alpha),
	                             Parameter = "alpha",
	                             Samples = object$pars$alpha), 
			dplyr::tibble(iteration_ix = 1:length(object$pars$sigma),
			              Parameter = "sigma",
			              Samples = object$pars$sigma)) 



	P <- ncol(object$spec$X[[1]])
	probs <- object$pars$pi
	num_penalties <- length(object$spec$S)
	clnms <- Reduce(c,
	                lapply(1:K,function(x) paste0("K: " , x," ",colnames(object$spec$X[[1]]) ))
	                )

	beta <- object$pars$beta
	colnames(beta) <- c(colnames(object$mf$X),
						clnms)
	delta <- beta[,colnames(object$mf$X)]

	beta <- beta[,clnms]
	clnm_k <- lapply(1:K,function(x) grep(x=clnms,pattern = paste0("K: ",x),value=T))

	beta <- lapply(clnm_k,function(k) beta[,k])
	beta <- abind::abind(beta,along=3)


	scales <- object$pars$scales
	colnames(scales) <- paste0("tau_",1:(K*num_penalties))
	colnames(probs) <- paste0("pi","_",1:K)

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
				spec = object$spec
				)
	if(!is.null(object$pars$subj_b)){
		## Do better post processing here
		out$subj_b <- object$pars$subj_b
		out$subj_D <- object$pars$subj_D
	}

    structure(out, class = c("stapDP"))
}

