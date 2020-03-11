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

#' plots pairwise probability clustering plot
#'
#' @export
#' @param x stapDP object
#' @param sort boolean asking whether sorting algorithm should be used to sort
#' @return ggplot plot object
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn
#' @importFrom stringr str_replace
#'
plot_pairs.stapDP <- function(x,sort = FALSE){

	### To pass R CMD Check
	Group_2 <- Group_1 <- Probability <- NULL
	###

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
