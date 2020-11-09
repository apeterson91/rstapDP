
#' Create a stapDP object
#'
#' @param object A list provided by the fdp_staplm(er) functions
#' @return A stapDP object
#'
stapDP <- function(object){

	 Samples <- Parameter <- Lower <- Upper <- medp <- iteration_ix <- 
		. <- Distance <- Median <- P <-  id <- NULL
	
	if(!is.null(object$pars$subj_b)){
		glmod <- object$mf$glmod
		b <- reformat_b(object$pars$subj_b,glmod)
		D <- reformat_D(object$pars$subj_D,glmod)
		mer <- TRUE
	}else{
		mer <- FALSE
	}


	K <- object$K
	spec <- object$spec

	probs <- object$pars$pi
	meds <- apply(probs,2,median)
	ix <- order(meds,decreasing=TRUE)
	probs <- probs[,ix]
	colnames(probs) <- paste0("K: ",1:K)
	create_ixmat_vec <- function(K,K_product,ix){
		as.numeric(Reduce(cbind,t(t(matrix(1:(K*K_product), nrow = K_product, ncol = K ))[ix,])  ))
	}

	if(has_bw(spec)){
		num_penalties <- length(spec$S[[1]])
		ixmat_scale <- create_ixmat_vec(K,num_penalties,ix)
		tau_b <- object$pars$tau_b[,ixmat_scale]
		tau_w <- object$pars$tau_w[,ixmat_scale]
		colnames(tau_b) <- paste0("tau_b_",1:(K*num_penalties))
		colnames(tau_w) <- paste0("tau_w_",1:(K*num_penalties))
		scales <- cbind(tau_b,tau_w)

	}else{
		num_penalties <- length(object$spec$S)
		ixmat_scale <- create_ixmat_vec(K,num_penalties,ix)
		scales <- object$pars$scales[,ixmat_scale]
		colnames(scales) <- paste0("tau_",1:(K*num_penalties))
	}


	nms <- Reduce(c,lapply(spec$X,colnames))
	clnms_k <- lapply(1:K,function(x) paste0("K: " , x," ",nms ))
	clnms <- Reduce(c,clnms_k)

	ixmat_coef <- create_ixmat_vec(K,ncol(spec$X[[1]]),ix)
	beta <- cbind(object$pars$beta[,1:ncol(spec$mf$X)],object$pars$beta[,(ncol(spec$mf$X)+1):(ncol(object$pars$beta))][,ixmat_coef])
	colnames(beta) <- c(colnames(spec$mf$X),
						clnms)
	delta <- beta[,colnames(spec$mf$X)]
	betamat <- beta[,clnms]
	beta <- lapply(clnms_k,function(k) betamat[,k])
	beta <- abind::abind(beta,along=3)
	dimnames(beta)[[2]] <- nms
	dimnames(beta)[[3]] <- paste0("K: ",1:K)

	alpha <- matrix(object$pars$alpha,ncol=1,nrow=length(object$pars$alpha))
	colnames(alpha) <- "alpha"
	sigma <- matrix(object$pars$sigma,ncol=1,nrow=nrow(alpha))
	colnames(sigma) <- "sigma"


	parmat <- cbind(delta,betamat,alpha,sigma,probs,scales)
	if(mer){
		parmat <- cbind(parmat,b,D)
	}
	parameter_summary <- get_summary(parmat)


	ys <- object$pars$yhat
	gd <- expand.grid(id =paste("V_",1:ncol(object$pars$yhat)),
					  iteration_ix = 1:length(object$pars$alpha))

	colnames(ys) <- paste("V_",1:ncol(object$pars$yhat))

	ys <- suppressWarnings(dplyr::as_tibble(ys,quiet=T)) %>% 
	  dplyr::mutate(iteration_ix = 1:dplyr::n()) %>% 
		tidyr::pivot_longer(dplyr::contains("V_"),names_to="id",values_to="Samples") %>% 
		dplyr::mutate(id = as.integer(stringr::str_replace(id,"V_","")))

	yhat <- dplyr::tibble(iteration_ix = as.integer(gd$iteration_ix), 
						  Parameter = "yrep",
						  id = as.integer(gd$id))

	yhat <- suppressMessages(yhat %>% dplyr::right_join(ys))
	yhat <- rbind(yhat,dplyr::tibble(iteration_ix = rep(0,length(spec$mf$y)),Parameter=rep("yobs",length(spec$mf$y)), 
									 id = 1:length(spec$mf$y),Samples = spec$mf$y))


	out <- list(beta = beta,
				summary = parameter_summary,
				delta = delta,
				alpha = alpha,
				sigma = sigma,
				probs = probs,
				yhat = yhat,
				scales = scales,
				pmat = object$pars$pmat,
				cmat = object$pars$clabels,
				model = list(formula = object$formula,
							 K=object$K,
							 y=spec$mf$y,
							 alpha_a = object$alpha_a,
							 alpha_b = object$alpha_b),
				spec = object$spec)

	if(!is.null(object$mf$glmod)){
		out$glmod <- glmod
		out$subj_b <- b
		out$subj_D <- D
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


get_summary <- function(parmat){
	nms <- colnames(parmat)
	mean <- colMeans(parmat)
	sd <- apply(parmat,2,sd)
	qs <- t(apply(parmat,2,function(x) quantile(x,c(0.1,.25,.5,.75,.9),na.rm=T)))
	n_eff <- apply(parmat,2,rstan::ess_tail)
	Rhat <- apply(parmat,2,rstan::Rhat)
	out <- cbind(mean,sd,qs,n_eff,Rhat)
}
