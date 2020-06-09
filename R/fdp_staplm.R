
#' Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @param formula Similar as for \code{\link[stats]{lm}}. 
#' @param knots quantiles of distance
#' @param data frame
#' @param tau_0 normal base measure scale
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param K truncation number
#' @param penalize boolean that denotes whether parameters are penalized or not
#' @param iter_max maximum number of iterations
#' @param burn_in number of burn in iterations
#' @param thin number by which to thin samples
#' @param seed rng initializer
#' 
#' @export
#' 
fdp_staplm <- function(formula,
					   knots,
					   subject_data,
					   dt_data = NULL,
					   w = rep(1,nrow(subject_data)),
					   tau_df = 1,
					   alpha_a,
					   alpha_b, 
					   K,
					   penalize = T,
					   iter_max,
					   burn_in,
					   thin = 1,
					   seed = NULL){

	stopifnot(burn_in<iter_max && burn_in > 0)
	stopifnot(alpha_a>0)
	stopifnot(alpha_b>0)
	stopifnot(thin>0)
	call <- match.call(expand.dots = TRUE)
	mf <- match.call(expand.dots = FALSE)
	mf$formula <- formula
	m <- match(c("formula"),table = names(mf), nomatch=0L)
	mf <- mf[c(1L,m)]
	mf$data <- subject_data
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf,parent.frame())
	mt <- attr(mf,"terms")
	if(is.empty.model(mt))
		stop("No intercept or predictors specified.",.call = FALSE)


	y <- model.response(mf, "numeric")
	Z <- model.matrix(mt,mf)
	X <- splines::bs(subject_data$Distances/max(subject_data$Distances))
	S <- diag(ncol(Z) + ncol(X)*K)
	S[1:ncol(Z),1:ncol(Z)] <- 0 ## improper priors on Z beta's
	# X <- subject_data %>% dplyr::left_join(dt_data) %>% split(.$id) %>% 
	#   purrr::map(.,function(x){
	#     if(all(is.na(x$Distance)))
	#       return(matrix(0,ncol=4))
	#     else
	#       return(colSums(cbind(1,x$Distance/length(Distance))))
	#   }) %>% do.call(rbind,.)
	# D <- get_D(knots)
	# X <- get_X_transformed(X,D)

	
	fit <- fdp_staplm.fit(y,Z,X, S, alpha_a,alpha_b,K,penalize,tau_df,w,iter_max,burn_in,thin,seed)
    fit <- list(beta = fit$beta,
                probs = fit$pi,
                sigma = fit$sigma,
                alpha = fit$alpha,
				yhat = fit$yhat,
				cluster_mat = fit$cluster_assignment,
                scales = if(penalize) fit$tau else matrix(rep(tau_df,K),ncol=K),
                pmat = fit$PairwiseProbabilityMat,
                clabels = fit$cluster_assignment,
                Znames = colnames(Z),
                ncol_X = ncol(X),
                penalize = penalize,
                formula = formula,
				y = y,
                K = K)
	return(stapDP(fit))
}


#'  Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @param y vector of outcomes
#' @param Z design matrix
#' @param X stap design matrix
#' @param S desi
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param K truncation number for DP mixture components
#' @param nu_0 inverse-chi-square degree of freedom for tau base measure 
#' @param w weights for weighted regression - default is 1
#' @param iter_max maximum number of iterations
#' @param burn_in number of burn in iterations
#' @param thin number by which to thin samples
#' @param seed rng initializer
#' @export
#' 
fdp_staplm.fit <- function(y,Z,X,S,
						   alpha_a,
						   alpha_b, 
						   K,
						   penalize,
						   tau_df,
						   w = rep(1,length(y)),
						   iter_max,burn_in,
						   thin,seed = NULL){

	stopifnot(tau_df>0)
	stopifnot(ncol(S) == nrow(S))
	stopifnot(ncol(S) == ncol(Z) + ncol(X)*K)
  if(is.null(seed)){
    seed <- 3413
  }
  num_posterior_samples <- length(seq(from=burn_in+1,to=iter_max,by=thin))

  if(penalize)
	  fit <- stappDP_fit(y,Z,X,S,w,tau_df,alpha_a,alpha_b,K,iter_max,burn_in,thin,seed,num_posterior_samples)
  else
	  fit <- stapDP_fit(y,Z,X,w,
						  tau_df,alpha_a,alpha_b,
						  K,iter_max,burn_in,thin,
						  seed,num_posterior_samples)


  return(fit)
}
