
#' Create a stapDP object
#'
#' @param object A list provided by the fdp_staplm.fit function
#' @return A stapDP object
#'
stapDP <- function(object){

	K <- Samples <- Parameter <- Lower <- Upper <- medp <- iteration_ix <- 
		. <- Distance <- Median <- P <-  id <- NULL
	pardf <- rbind(dplyr::tibble(iteration_ix = 1:length(object$alpha),
	                             Parameter = "alpha",
	                             Samples = object$alpha), 
			dplyr::tibble(iteration_ix = 1:length(object$sigma),
			              Parameter = "sigma",
			              Samples = object$sigma)) 

	fixef <- object$beta[,1:length(object$Znames)]

	ranef <- object$beta[,(length(object$Znames)+1):ncol(object$beta)]

	colnames(fixef) <- object$Znames

	gd <- expand.grid(p=1:object$ncol_X,k=1:object$K)

	colnames(ranef) <- paste0("Distance K:",gd$k,"_P:",gd$p)

	fixef <- dplyr::as_tibble(fixef,quiet=T) %>% 
	  dplyr::mutate(iteration_ix = 1:dplyr::n()) %>% 
	  tidyr::gather(object$Znames,key = "Parameter", value = "Samples")

	ranef <- dplyr::as_tibble(ranef,quiet=T) %>% dplyr::mutate(iteration_ix = 1:dplyr::n()) %>% 
	  tidyr::gather(colnames(ranef),key="Parameter",value="Samples") %>% 
	  dplyr::mutate(K = as.integer(stringr::str_replace(stringr::str_extract(Parameter,"K:[0-9]{2}|K:[0-9]{1}"),"K:","")),
					P = stringr::str_replace(stringr::str_extract(Parameter,"P:[0-9]{2}|P:[0-9]{1}"),"P:",""))

	scales <- object$scales
	colnames(scales) <- paste0("tau_",1:(object$K))
	probs <- dplyr::as_tibble(object$probs,quiet=T) 
	colnames(probs) <- paste0("pi","_",1:object$K)
	probs <- probs %>% dplyr::mutate(iteration_ix = 1:dplyr::n()) %>% 
		tidyr::gather(dplyr::contains("pi"),key="Parameter",value="Samples") %>% 
		dplyr::mutate(K = as.integer(stringr::str_replace( stringr::str_extract(Parameter,"_[0-9]{2}|_[0-9]{1}"),"_","")))

	ys <- object$yhat
	gd <- expand.grid(id =paste("V_",1:ncol(object$yhat)),
					  iteration_ix = 1:length(object$alpha))

	colnames(ys) <- paste("V_",1:ncol(object$yhat))

	ys <- dplyr::as_tibble(ys,quiet=T) %>% dplyr::mutate(iteration_ix = 1:dplyr::n()) %>% 
		tidyr::gather(dplyr::contains("V_"),key="id",value="Samples") %>% 
		dplyr::mutate(id = as.integer(stringr::str_replace(id,"V_","")))

	yhat <- dplyr::tibble(iteration_ix = as.integer(gd$iteration_ix), 
						  Parameter = "yrep",
						  id = as.integer(gd$id))

	yhat <- suppressMessages(yhat %>% dplyr::right_join(ys))
	yhat <- rbind(yhat,dplyr::tibble(iteration_ix = rep(0,length(object$y)),Parameter=rep("yobs",length(object$y)), id = 1:length(object$y),Samples = object$y))


	out <- list(pardf = pardf,
				fixef = fixef,
				ranef = ranef,
				probs = probs,
				yhat = yhat,
				scales = scales,
				pmat = object$pmat,
				cmat = object$cluster_mat,
				model = list(formula = object$formula,
							 K=(object$K),
							 y=object$y,
							 sobj = object$sobj)
				)

    structure(out, class = c("stapDP"))
}

