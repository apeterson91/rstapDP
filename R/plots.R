#' Plotting functions for stapDP objects 

#' Plots pairwise probability clustering plot
#' 
#' @template rodriguez
#'
#' @export
#' @param x stapDP object
#' @param sort boolean asking whether sorting algorithm should be used to sort
#' pairwise probablity
#' @param sample number of subjects to sample in order to speed computation
#' @return ggplot plot object
#' @seealso the supplementary section of the reference for the sorting algorithm.
plot_pairs <- function(x,sort = FALSE, sample = 0)
    UseMethod("plot_pairs")

#'
#'
#' @export
#' @describeIn plot_pairs plot_pairwise probability matrix
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#' @importFrom stringr str_replace
#'
plot_pairs.stapDP <- function(x,sort = FALSE,sample = 0){

	### To pass R CMD Check
	Group_2 <- Group_1 <- Probability <- NULL
	###

	if(sample>0){
		ics <- sample(1:nrow(x$pmat),sample,replace=F)
		P <- x$pmat[ics,ics]
	}else
		P <- x$pmat

	makeSymm <- function(m) {
	  m[upper.tri(m)] <- t(m)[upper.tri(m)]
	  return(m)
	}

	P <- makeSymm(P)
	if(sort){

		i <- which(P == max(P),arr.ind = T)[1]
		j <- which(P == max(P),arr.ind = T)[2]

		A <- c(i,j)
		Omega <- 1:nrow(P)
		A_c <- setdiff(Omega,A)
		while(length(A_c)){
		  probs <- sapply(A_c,function(x) max(P[x,A]))
		  j <- A_c[which.max(probs)]
		  A <- c(A,j)
		  A_c <- setdiff(Omega,A)
		}
	}
	else{
		A <- 1:nrow(P)
	}


	p <- suppressMessages(dplyr::as_tibble(P[A,A])) %>% dplyr::mutate(Group_1 = 1:dplyr::n()) %>%
	  tidyr::gather(dplyr::contains("V"), key = "Group_2", value = "Probability") %>%
	  dplyr::mutate("Group_2" = as.numeric(stringr::str_replace(Group_2,"V",""))) %>%
	  ggplot(aes(x=Group_1,y=Group_2,fill=Probability)) +
	  geom_tile() + scale_fill_gradientn(colours=c("white","grey","black"),limits=c(0,1)) +
	  ggplot2::theme_bw() + ggplot2::labs(title = "Pairwise Probability of Function Clustering",
										  x="Subject 1", y = "Subject 2") +
	  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
			panel.grid.minor = ggplot2::element_blank())


    return(p)
}



#' Cluster-Specific Spatial Temporal Effects
#' 
#' @export
#' @param x stapDP object
#' @param p probability contained in credible interval
#' @param style one of "color" or "facet" for different plotting options
#' @param  prob_filter all mixture components with median probability < prob_filter are excluded from the plot
#' @param ... ignored
#' @return plot with cluster effect across space
#' 
plot.stapDP <- function(x,p = 0.95, 
						style = "color",
						prob_filter = 0.1,
						...){

	K <- Samples <- Parameter <- Lower <- Upper <- medp <- iteration_ix <- 
		Model <- Prob <- . <- Distance <- Median <- P <- y <- RSS <- mnRSS <- NULL

	stopifnot(p>=0 && p<=1)
	l <-  .5 - p/2
	u <- .5 + p/2
	spec <- x$spec
	term <- x$spec$term[1]
	comp <- x$spec$component[1]
	gd_eta <- get_stap(spec,term,comp,x$beta)

	kprob <- apply(x$probs,2,median)
	ks_to_keep <- which(kprob>prob_filter)

	if(has_bw(spec)){
		pltdf <- purrr::map_dfr(ks_to_keep,function(x)
								  dplyr::tibble(Distance = gd_eta$grid$Distance,
												K = x,
												Model = "Between",
												Lower = apply(gd_eta$eta_b[,,x],1,function(y) quantile(y,l)),
												Median = apply(gd_eta$eta_b[,,x],1,median),
												Upper = apply(gd_eta$eta_b[,,x],1,function(y) quantile(y,u))) %>%
								  rbind(.,
										dplyr::tibble(Distance = gd_eta$grid$Distance,
													  K = x,
													  Model = "Within",
													  Lower = apply(gd_eta$eta_w[,,x],1,function(y) quantile(y,l)),
													  Median = apply(gd_eta$eta_w[,,x],1,median),
													  Upper = apply(gd_eta$eta_w[,,x],1,function(y) quantile(y,u)))
								  )
		) 
		mid <- median(gd_eta$grid$Distance)
		max_eta <- quantile(gd_eta$eta_b[,,ks_to_keep],0.99)
		kprob <- dplyr::tibble(x = mid, 
		                       y = max_eta,
		                       K = factor(ks_to_keep),
		                       Model = factor("Between"),
		                       Prob = paste("P(.) = ",round(100*kprob[ks_to_keep],2)))
		p <- pltdf %>% dplyr::mutate(K = factor(K),
									 Model=factor(Model)) %>%
		ggplot2::ggplot(ggplot2::aes(x=Distance,y=Median,linetype=K)) + 
		  ggplot2::geom_line() +
		  ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower,ymax=Upper),alpha=0.3)+
		ggplot2::geom_hline(ggplot2::aes(yintercept=0),linetype=2,color='red') + 
		ggplot2::facet_wrap(~Model + K,nrow = 2,ncol=length(ks_to_keep)) +  
		  ggplot2::geom_label(ggplot2::aes(x=x,y=y,label=Prob),data=kprob) + 
		  ggplot2::theme(strip.background = ggplot2::element_blank()) + 
		  ggplot2::labs(y="Exposure Effect")
	return(p)
	}
	gd <- gd_eta$grid

	eta <- gd_eta$eta


	pltdf <- purrr::map_dfr(ks_to_keep,function(x)
							  dplyr::tibble(Distance = gd$Distance,
							                K = x,
							                Lower = apply(eta[,,x],1,function(y) quantile(y,l)),
							                Median = apply(eta[,,x],1,median),
							                Upper = apply(eta[,,x],1,function(y) quantile(y,u))))


	pltdf %>% dplyr::mutate(K = factor(K)) %>% 
	  ggplot2::ggplot(ggplot2::aes(x=Distance,y=Median,linetype=K)) + 
	  ggplot2::geom_hline(ggplot2::aes(yintercept = 0),linetype=2,color='red')-> p   

	if (style=="color"){
		p + ggplot2::geom_line(ggplot2::aes(color=K)) + ggplot2::theme_bw() + 
		  ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
		  ggplot2::labs(y="Exposure Effect") -> pl
	}else{
	  mid <- median(gd$Distance)
	  max_eta <- quantile(eta[,,ks_to_keep],0.99)
	  kprob <- dplyr::tibble(x = mid, 
	                         y = max_eta,
	                         K = factor(ks_to_keep),
	                         Prob = paste("P(.) = ",round(100*kprob[ks_to_keep],2)))
	  
		p + ggplot2::geom_line() + ggplot2::theme_bw() + 
		  ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
		  ggplot2::facet_wrap(~K) + 
		  ggplot2::theme(strip.background=ggplot2::element_blank()) +
		  ggplot2::geom_label(aes(x=x,y=y,label=Prob),data=kprob) +
		  ggplot2::labs(y="Exposure Effect") -> pl
	}

	return(pl)

}

#' Diagnostic Traceplots
#' 
#' @export
#' @param x a stapDP object
#' @param par string of parameter to use
#' default is .1
#' 
traceplots <- function(x,par = NULL){
	UseMethod("traceplots")
}

#'
#' 
#' @export
#' @describeIn traceplots
#' 
traceplots.stapDP <- function(x,par = NULL){

	mat <- as.matrix(x)
	if(!is.null(par)){
		stopifnot(par %in% colnames(mat))
	}
	
	p <- bayesplot::mcmc_trace(list(as.matrix(x)[,par,drop=F]))
	return(p)
}


#' Posterior Predictive Checks 
#'
#' @export
#' @template reference-bda
#' @param x a stapDP object
#' @param num_reps number of yhat samples to plot
#'
ppc <- function(x,num_reps = 20)
	UseMethod("ppc")

#' Posterior Predictive Checks
#'
#' @export
#' @describeIn ppc
#'
ppc.stapDP <- function(x,num_reps = 20){

	Samples <- Parameter <- iteration_ix <- NULL

	max_num_iter <- max(x$yhat$iteration_ix)
	if(num_reps > max_num_iter)
		stop("Not enough iterations to display this many samples")

	samp <- c(0,sample(1:max(x$yhat$iteration_ix),num_reps))
	p <- suppressWarnings(x$yhat %>% dplyr::filter(iteration_ix %in% samp) %>% 
		ggplot2::ggplot(ggplot2::aes(x=Samples,color=Parameter,group=iteration_ix,alpha=Parameter)) + 
		ggplot2::geom_density() + ggplot2::theme_bw()+ ggplot2::theme(legend.title=ggplot2::element_blank()) +   
		ggplot2::scale_colour_manual(values=c("black","grey")) + ggplot2::scale_alpha_discrete(range=c(1,.1)))

	return(p)
}
