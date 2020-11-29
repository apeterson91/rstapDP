
#' Create a stapDP object
#'
#' @param object A list provided by the fdp_staplm(er) functions
#' @return A stapDP object
#' @importFrom utils tail
#'
stapDP <- function(object){

	 Samples <- Parameter <- Lower <- Upper <- medp <- iteration_ix <- 
		. <- Distance <- Median <- P <-  id <- NULL
	

	## for sorting matrices by cluster probability
	create_ixmat_vec <- function(K,K_product,ix){
		as.numeric(Reduce(cbind,t(t(matrix(1:(K*K_product), nrow = K_product, ncol = K ))[ix,])  ))
	}

	collapse_pars <- function(par_name,K,K_product,ics){
		ixmats <- lapply(ics,function(x) create_ixmat_vec(K,K_product,x))
		par <- Reduce(rbind,purrr::map2(pars,ixmats,function(x,y) x[[par_name]][,y]))
	}


	K <- object$K
	spec <- object$spec
	pars <- object$pars

	if(!is.null(pars[[1]]$subj_b)){
		glmod <- object$mf$glmod
		b <- Reduce(rbind,lapply(pars,function(x) reformat_b(x$subj_b,glmod)))
		D <- Reduce(rbind,lapply(pars, function(x) reformat_D(x$subj_D,glmod)))
		mer <- TRUE
	}else{
		mer <- FALSE
	}

	ix <- lapply(pars,function(x) order(apply(x$pi,2,median),decreasing=TRUE))
	probs <- Reduce(rbind,purrr::map2(pars,ix,function(x,y) x$pi[,y]))
	colnames(probs) <- paste0("pi_K: ",1:K)

	if(has_bw(spec)){
		num_penalties <- length(spec$S[[1]])
		tau_b <- collapse_pars("tau_b",K,num_penalties,ix)
		tau_w <- collapse_pars("tau_w",K,num_penalties,ix)
		colnames(tau_b) <- paste0("tau_b_",1:(K*num_penalties))
		colnames(tau_w) <- paste0("tau_w_",1:(K*num_penalties))
		scales <- cbind(tau_b,tau_w)
	}else{
		num_penalties <- length(object$spec$S)
		scales <- collapse_pars("scales",K,num_penalties,ix)
		colnames(scales) <- paste0("tau_",1:(K*num_penalties))
	}


	nms <- Reduce(c,lapply(spec$X,colnames))
	clnms_k <- lapply(1:K,function(x) paste0("K: " , x," ",nms ))
	clnms <- Reduce(c,clnms_k)
	delta_ics <- 1:ncol(spec$mf$X)
	beta_ics <- (tail(delta_ics,1)+1):ncol(pars[[1]]$beta)
	beta_prod <- ncol(spec$X[[1]]) + has_bw(spec)*ncol(spec$X[[1]])

	bixmats <- lapply(ix,function(x) create_ixmat_vec(K,beta_prod,x))
	betamat <- Reduce(rbind,purrr::map2(pars,bixmats,function(x,y) x$beta[,beta_ics][,y]))
	delta <- Reduce(rbind,lapply(pars,function(x) x$beta[,delta_ics]))
	colnames(delta) <- colnames(spec$mf$X)
	colnames(betamat) <- clnms
	beta <- lapply(clnms_k,function(k) betamat[,k])
	beta <- abind::abind(beta,along=3)
	dimnames(beta)[[2]] <- nms
	dimnames(beta)[[3]] <- paste0("K: ",1:K)

	alpha <- Reduce(rbind,lapply(pars, function(x) matrix(x$alpha,ncol=1,nrow=length(x$alpha))))
	colnames(alpha) <- "alpha"
	sigma <- Reduce(rbind,lapply(pars, function(x) matrix(x$sigma,ncol=1,nrow=length(x$sigma))))
	colnames(sigma) <- "sigma"


	parmat <- cbind(delta,betamat,alpha,sigma,probs,scales)
	if(mer){
		parmat <- cbind(parmat,b,D)
	}
	parameter_summary <- get_summary(parmat,length(pars))


	ys <- Reduce(rbind,lapply(pars, function(x) x$yhat))
	gd <- expand.grid(id =paste0("V_",1:ncol(pars[[1]]$yhat)),
					  iteration_ix = 1:nrow(alpha))

	colnames(ys) <- paste("V_",1:ncol(pars[[1]]$yhat))

	ys <- dplyr::as_tibble(ys,quiet=T) %>% 
	  dplyr::mutate(iteration_ix = 1:dplyr::n()) %>% 
		tidyr::pivot_longer(dplyr::contains("V_"),names_to="id",values_to="Samples") %>% 
		dplyr::mutate(id = as.integer(stringr::str_replace(id,"V_","")))

	yhat <- dplyr::tibble(iteration_ix = as.integer(gd$iteration_ix), 
						  Parameter = "yrep",
						  id = as.integer(gd$id))

	yhat <- yhat %>% dplyr::right_join(ys,by=c('iteration_ix','id'))
	yhat <- rbind(yhat,dplyr::tibble(iteration_ix = rep(0,length(spec$mf$y)),
									 Parameter=rep("yobs",length(spec$mf$y)), 
									 id = 1:length(spec$mf$y),Samples = spec$mf$y))


	out <- list(beta = beta,
				summary = parameter_summary,
				delta = delta,
				alpha = alpha,
				sigma = sigma,
				probs = probs,
				yhat = yhat,
				scales = scales,
				pmat = Reduce(`+`,lapply(pars,function(x) x$pmat))/length(pars),
				cmat = Reduce(rbind,lapply(pars,function(x) x$clabels)),
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
create_chain_mat <- function(mat_col,num_samples){
	ixmat <- cbind(seq(from = 1L, to = length(mat_col),by=num_samples),
				   seq(from = num_samples,to=length(mat_col),by=num_samples))
	mat <- Reduce(cbind,lapply(1:nrow(ixmat),function(x) mat_col[ixmat[1,1]:ixmat[1,2]]))
}


get_summary <- function(parmat,chains){
	nms <- colnames(parmat)
	num_samples <- nrow(parmat) / chains
	mean <- colMeans(parmat)
	sd <- apply(parmat,2,sd)
	qs <- t(apply(parmat,2,function(x) quantile(x,c(0.1,.25,.5,.75,.9),na.rm=T)))
	cols <- lapply(1:ncol(parmat),function(x) create_chain_mat(parmat[,x],num_samples))
	n_eff <- Reduce(rbind,lapply(cols,rstan::ess_tail))
	Rhat <- Reduce(rbind,lapply(cols,rstan::Rhat)) 
	out <- cbind(mean,sd,qs,n_eff,Rhat)
	colnames(out)[8:9] <- c("n_eff","Rhat")
	return(out)
}
