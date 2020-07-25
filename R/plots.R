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


	p <- dplyr::as_tibble(P[A,A]) %>% dplyr::mutate(Group_1 = 1:dplyr::n()) %>%
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
#' @param switch one of "color" or "facet" for different plotting options
#' @param  prob_filter all mixture components with median probability < prob_filter are excluded from the plot
#' @param ... ignored
#' @return plot with cluster effect across space
#' 
plot.stapDP <- function(x,p = 0.95, 
						switch = "color",
						prob_filter = 0.1,...){

	K <- Samples <- Parameter <- Lower <- Upper <- medp <- iteration_ix <- 
		. <- Distance <- Median <- P <- y <- RSS <- mnRSS <- NULL
	gd <- data.frame(var = seq(from = 0, to = 1, by = 0.01))
	colnames(gd) <- x$model$sobj$term

	mat <- mgcv::Predict.matrix(x$model$sobj,gd)

	ks <- x$probs %>% dplyr::group_by(K) %>% dplyr::summarise(medp=median(Samples))%>%
		dplyr::filter(medp>prob_filter) %>% dplyr::pull(K)

	if(length(ks)==0)
		stop("No clusters included with set prob_filter, try a lower")
	P <- max(as.integer(x$ranef$P))

	x$ranef %>% dplyr::filter(K %in% ks) %>% tidyr::spread(P,Samples) %>% 
	  dplyr::mutate_if(is.double,function(x) tidyr::replace_na(x,0)) %>%
	  dplyr::group_by(iteration_ix,K) %>% 
	  dplyr::summarise_if(is.double,sum) %>% 
	  dplyr::ungroup() %>% 
	  dplyr::select(-iteration_ix) %>% 
	  split(.$K) %>% 
	  purrr::map2_dfr(.,names(.),function(x,y) {
		mt <- mat %*% (x %>% dplyr::select(as.character(1:P)) %>% as.matrix() %>% t())  
		colnames(mt) <- paste0("ix_",1:ncol(mt))
		
		df <- dplyr::as_tibble(mt) %>% dplyr::mutate(K= rep(y,dplyr::n()),
													 Distance = gd[,1]) %>% 
		  tidyr::gather(dplyr::contains("ix_"),key="Iteration_ix",value="Samples")
		}) -> pltdf
	pltdf %>% dplyr::group_by(Distance,K) %>% 
	  dplyr::summarise(Lower = quantile(Samples,0.025),
					   Median = median(Samples),
					   Upper = quantile(Samples,0.975)) %>% 
	  ggplot2::ggplot(ggplot2::aes(x=Distance,y=Median,linetype=K)) -> p   
	if (switch=="color"){
		p + ggplot2::geom_line(ggplot2::aes(color=K)) + ggplot2::theme_bw() + 
		  ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
		  ggplot2::labs(y="Exposure Effect") -> pl
	}else{
		p + ggplot2::geom_line() + ggplot2::theme_bw() + 
		  ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),alpha=0.3) + 
		  ggplot2::facet_wrap(~K) + 
		  ggplot2::labs(y="Exposure Effect") -> pl
	} 	

	return(pl)

}

#' Diagnostic Traceplots
#' 
#' @export
#' @param x a stapDP object
#' @param par string of parameter to use
#' @param prob_filter median probability of cluster components to include in plot
#' default is .1
#' 
traceplots <- function(x,par=c("probs"),prob_filter = .1){
	UseMethod("traceplots")
}

#'
#' 
#' @export
#' @describeIn traceplots
#' 
traceplots.stapDP <- function(x,par=c("probs"),prob_filter = .1){
	stopifnot(par %in% c("probs","fixef","ranef","alpha","sigma"))

	K <- Samples <- Parameter <- Lower <- Upper <- medp <- iteration_ix <- NULL

	ks <- x$probs %>% dplyr::group_by(K) %>% dplyr::summarise(medp=median(Samples))%>%
		dplyr::filter(medp>prob_filter) %>% dplyr::pull(K)

	if(par %in% c("alpha","sigma"))
		p <- x$pardf %>% dplyr::filter(Parameter == !!par) %>% ggplot2::ggplot(ggplot2::aes(x=iteration_ix,y=Samples,color=Parameter)) + ggplot2::geom_line() + ggplot2::theme_bw()
	else{
		if(par == "fixef")
			p <- x[[par]] %>% ggplot2::ggplot(ggplot2::aes(x = iteration_ix ,y = Samples,color = Parameter)) + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::facet_wrap(~Parameter, scales = "free") + ggplot2::theme(strip.background=ggplot2::element_blank()) else p <- x[[par]] %>% dplyr::filter(K %in% ks) %>% 
				ggplot2::ggplot(ggplot2::aes(x = iteration_ix,y = Samples,color=Parameter)) + 
				ggplot2::geom_line() + 
				ggplot2::theme_bw() + ggplot2::facet_wrap(~Parameter + K,scales="free") + ggplot2::theme(strip.background = ggplot2::element_blank())
	}
	return(p)
}

#' Parameter Histograms
#'
#' @export
#' @param x a stapDP object
#' @param par string of parameter to use
#' @param prob_filter median probability of cluster components to include in plot
#' default is .1
#' 
plotpars <- function(x,par=c("probs"),prob_filter = .1)
	UseMethod("plotpars")

#' Parameter Histograms
#'
#' @export
#' @describeIn plotpars
#'
plotpars.stapDP <- function(x,par=c("probs"),prob_filter = .1){

	stopifnot(par %in% c("probs","fixef","ranef","alpha","sigma"))

	ks <- x$probs %>% dplyr::group_by(K) %>% dplyr::summarise(medp=median(Samples))%>%
		dplyr::filter(medp>prob_filter) %>% dplyr::pull(K)

	K <- Samples <- Parameter <- Lower <- Upper <- medp <- NULL

	if(par %in% c("alpha","sigma"))
		p <- x$pardf %>% dplyr::filter(Parameter == !!par) %>% ggplot2::ggplot(ggplot2::aes(x=Samples,fill=Parameter)) + ggplot2::geom_histogram() + ggplot2::theme_bw()
	else{
		if(par=='fixef')
			p <- x[[par]] %>% ggplot2::ggplot(ggplot2::aes(x = Samples,fill=Parameter)) + ggplot2::geom_histogram() + 
				ggplot2::theme_bw() + ggplot2::facet_wrap(~Parameter,scales="free") + ggplot2::theme(strip.background = ggplot2::element_blank())
		else
			p <- x[[par]] %>% dplyr::filter(K %in% ks) %>%  ggplot2::ggplot(ggplot2::aes(x = Samples,fill=Parameter)) + ggplot2::geom_histogram() + 
				ggplot2::theme_bw() + ggplot2::facet_wrap(~Parameter + K,scales="free") + ggplot2::theme(strip.background = ggplot2::element_blank())
	}
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
