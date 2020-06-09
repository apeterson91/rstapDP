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

#' plot cluster effects plots
#' 
#' @export
#' @param x stapDP object
#' @param p probability contained in credible interval
#' @param switch one of "color" or "facet" for different plotting options
#' @param  prob_filter all mixture components with median probability < prob_filter are excluded from the plot
#' @return plot with cluster effect across space
#' 
plot_cluster_effects <- function(x, p = 0.95, switch = "color", prob_filter = 0.1)
    UseMethod("plot_cluster_effects")

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


#' plot cluster effects plots
#' 
#' @export
#' @describeIn plot_cluster_effects 
#' 
plot_cluster_effects.stapDP <- function(x,p = 0.95, switch = "color",prob_filter = 0.1){

	stop("To be Implimented")

}
#' Diagnostic Traceplots
#' 
#' @export
#' 
traceplots <- function(x,par=c("probs")){
	UseMethod("traceplots")
}


#' 
#' @export
#' 
traceplots.stapDP <- function(x,par=c("probs")){
	stopifnot(par %in% c("probs","fixef","ranef","alpha","sigma"))

	if(par %in% c("alpha","sigma"))
		p <- x$pardf %>% dplyr::filter(Parameter == !!par) %>% ggplot2::ggplot(ggplot2::aes(x=iteration_ix,y=Samples,color=Parameter)) + ggplot2::geom_line() + ggplot2::theme_bw()
	else{
		p <- x[[par]] %>% ggplot2::ggplot(ggplot2::aes(x = iteration_ix,y = Samples,color=Parameter)) + ggplot2::geom_line() + 
			ggplot2::theme_bw() + ggplot2::facet_wrap(~Parameter,scales="free") + ggplot2::theme(strip.background = ggplot2::element_blank())
	}
	return(p)
}

#' Parameter Histograms
#'
#' @export
#' 
plotpars <- function(x,par=c("probs"))
	UseMethod("plotpars")

#' Parameter Histograms
#'
#' @export
#'
plotpars.stapDP <- function(x,par=c("probs")){
	if(par %in% c("alpha","sigma"))
		p <- x$pardf %>% dplyr::filter(Parameter == !!par) %>% ggplot2::ggplot(ggplot2::aes(x=Samples,fill=Parameter)) + ggplot2::geom_histogram() + ggplot2::theme_bw()
	else{
		p <- x[[par]] %>% ggplot2::ggplot(ggplot2::aes(x = Samples,fill=Parameter)) + ggplot2::geom_histogram() + 
			ggplot2::theme_bw() + ggplot2::facet_wrap(~Parameter,scales="free") + ggplot2::theme(strip.background = ggplot2::element_blank())
	}
	return(p)

}

#' Posterior Predictive Checks 
#'
#'@export
#'
ppc <- function(x,num_reps = 20)
	UseMethod("ppc")

#' Posterior Predictive Checks
#'
#' @export
#'
ppc.stapDP <- function(x,num_reps = 20){

	samp <- c(0,sample(1:max(x$yhat$iteration_ix),num_reps))
	p <- suppressWarnings(x$yhat %>% dplyr::filter(iteration_ix %in% samp) %>% 
		ggplot2::ggplot(ggplot2::aes(x=Samples,color=Parameter,group=iteration_ix,alpha=Parameter)) + 
		ggplot2::geom_density() + ggplot2::theme_bw()+ ggplot2::theme(legend.title=ggplot2::element_blank()) +   
		ggplot2::scale_colour_manual(values=c("black","grey")) + ggplot2::scale_alpha_discrete(range=c(1,.1)))

	return(p)
}
