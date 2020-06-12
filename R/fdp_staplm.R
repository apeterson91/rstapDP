
#' Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @param formula Similar as for \code{\link[stats]{lm}}. 
#' @param subject_data data frame of subject level measurements
#' @param stap_term call to \code{s_stap} with term in dt_data
#' @param dt_data data frame containing distances and/or times associated with subject data
#' @param w weights for weighted regression - default is vector of ones 
#' @param tau_df if penalize = F: normal base measure scale 
#'        if penalize true this is the degrees of freedom for the inverse chi square base measure. 
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param K truncation number
#' @param penalize boolean that denotes whether parameters are penalized or not
#' @param iter_max maximum number of iterations
#' @param burn_in number of burn in iterations
#' @param thin number by which to thin samples
#' @param knots optional nots argument passed to \code{\link[mgcv]{smooth.construct}}.
#' @param id_colname string of id column name that is the unique id in both subject_data and dt_data
#' @param seed random number generator seed will be set to default value if not by user
#' 
#' @importFrom stats is.empty.model model.matrix model.response
#' @export
#' 
fdp_staplm <- function(formula,
					   subject_data,
					   stap_term,
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
					   knots = NULL,
					   id_colname = 'id',
					   seed = NULL){

	## Parameter check
	stopifnot(burn_in<iter_max && burn_in > 0)
	stopifnot(alpha_a>0)
	stopifnot(alpha_b>0)
	stopifnot(thin>0)
	stopifnot(id_colname %in% colnames(dt_data))
	stopifnot(id_colname %in% colnames(subject_data))
	stopifnot(stap_term$term %in% colnames(dt_data))
	stopifnot(length(id_colname)==1)
	## 

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
	## Handle Zero exposure in dt_data
	jdf <- dplyr::left_join(subject_data,dt_data,by=id_colname) %>% 
	  dplyr::mutate_if(is.double,function(x) tidyr::replace_na(x,0))
	foo <- mgcv::smooth.construct(stap_term,data=jdf,knots=knots)
	M <- Matrix::fac2sparse(as.factor(jdf[,id_colname,drop=T]))
	X <- as.matrix(M %*% foo$X)
	S <- kronecker(diag(K),foo$S[[1]])
	S <- rbind(matrix(0,ncol=ncol(S)+ncol(Z),nrow=ncol(Z)),
	           cbind(matrix(0,nrow=nrow(S),ncol=ncol(Z)),S))
	
	
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
                sobj = foo,
                y = y,
                K = K)
	return(stapDP(fit))
}


#'  Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @param y vector of outcomes
#' @param Z design matrix
#' @param X stap design matrix
#' @param S penalty matrix from (e.g.) \code{\link[mgcv]{smooth.construct}} 
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param K truncation number for DP mixture components
#' @param penalize boolean value indicating whether or not the stap parameters
#' should be penalized via a normal prior/L2 penalty (as opposed to using an improper prior)
#' @param tau_df if penalize = F: normal base measure scale 
#'        if penalize true this is the degrees of freedom for the inverse chi square base measure. 
#' @param w weights for weighted regression - default is vector of ones 
#' @param iter_max maximum number of iterations
#' @param burn_in number of iterations to burn-in
#' @param thin number by which to thin samples
#' @param seed random number generator seed will be set to default value if not by user
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
	stopifnot(length(w) == length(y))
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
