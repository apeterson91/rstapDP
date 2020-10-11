#' Plotting functions for stapDP objects 


#'
#' @importFrom bendr plot_pairs
#' @export
plot_pairs.stapDP <- function(x,sample = NULL, sort = FALSE)
	UseNextMethod("plot_pairs")


#' @export
bendr::plot_pairs




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
	
	p <- bayesplot::mcmc_trace(list(mat[,par,drop=F]))
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
