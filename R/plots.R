#' plots pairwise probability clustering plot
#' 
#' @template rodriguez
#'
#' @export
#' @param x stapDP object
#' @param sort boolean asking whether sorting algorithm should be used to sort
#' pairwise probablity
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#' @return ggplot plot object
#' @seealso the supplementary section of the reference for the sorting algorithm.
#'

plot_pairs <- function(x,sort = FALSE)
    UseMethod("plot_pairs")

#' plot cluster effects plots
#' 
#' @export
#' @param x stapDP object
#' 
plot_cluster_effects <- function(x, p = 0.95, switch = "color", prob_filter = 0.1)
    UseMethod("plot_cluster_effects")

#' plots pairwise probability clustering plot
#'
#' @export
#' @param x stapDP object
#' @param sort boolean asking whether sorting algorithm should be used to sort
#' @return ggplot plot object
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
#' @param x stapDP object
#' 
plot_cluster_effects.stapDP <- function(x,p = 0.95, switch = "color",prob_filter = 0.1){

	if(!(switch %in% c("color","facet"))){
		stop("switch must be one of `color` or `facet` ")
	}
	stopifnot(p>0 & p<1)
	K <- ncol(x$pi)
	P_two <- length(x$X_ranges)
	start <- seq(from = 3, to = (2+(K*P_two)), by = P_two)
	end <- start+(P_two-1)
	grids <- seq(from = x$X_ranges[[2]][1],to=x$X_ranges[[2]][2], length.out = 150)
	grids <- cbind(1,grids,grids^2)
	lo <- (1 - p)/2
	hi <-  p + lo
	los <- lapply(1:K, function(k) apply(tcrossprod(grids, x$beta[,start[k]:end[k],drop=F]),1,function(x) quantile(x,lo)))
	his <- lapply(1:K, function(k) apply(tcrossprod(grids, x$beta[,start[k]:end[k],drop=F]),1,function(x) quantile(x,hi)))
	meds <- lapply(1:K, function(k) apply(tcrossprod(grids, x$beta[,start[k]:end[k],drop=F]),1,function(x) quantile(x,.5)))

	out <- purrr::map_dfr(1:K,function(k) {
	  dplyr::tibble(K = k,
	                Grid = grids[,2],
	                med = meds[[k]],
	                lower = los[[k]],
	                upper = his[[k]])
	})
	pi_df <- dplyr::tibble(K = 1:K,
	                pi = apply(fit$pi,2,median))

	if(switch=="color")
		p <- out %>% dplyr::left_join(pi_df) %>% dplyr::filter(pi>prob_filter) %>% 
	  dplyr::mutate(K= factor(K)) %>% 
			ggplot2::ggplot(aes(x=Grid,y=med,fill=K)) + 
			ggplot2::geom_line() + ggplot2::theme_bw() + 
			ggplot2::geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3)
	else
		p <- out %>% dplyr::left_join(pi_df) %>% dplyr::filter(pi>prob_filter) %>% 
	  dplyr::mutate(K= factor(K)) %>% 
			ggplot2::ggplot(aes(x=Grid,y=med)) + 
			ggplot2::geom_line() + ggplot2::theme_bw() + 
			ggplot2::geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3) + 
			ggplot2::facet_wrap(~K)

	p <- p + ggplot2::ggtitle("Cluster Spatial Effects") + ggplot2::geom_hline(aes(yintercept=0),linetype=2,color='red') +
	  ggplot2::xlab("Grid") + ggplot2::ylab("Change in g(E[Y])")
	return(p)

}
