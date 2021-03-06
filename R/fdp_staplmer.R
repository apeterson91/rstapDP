#' Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Mixed Effects Regression Model
#' 
#' @details This function fits a linear mixed effects regression model in a Bayesian paradigm with
#' improper priors assigned to the "standard" regression covariates designated 
#' in the formula argument and a Dirichlet process prior with normal-gamma base measure 
#' assigned to the stap basis function expansion using penalized splines via \code{\link[mgcv]{jagam}}.
#' normal priors are placed on the latent group variables and an improper prior is placed on the 
#' group covariance matrix leading to a Wishart posterior.
#'
#' The concentration parameter is assigned a gamma prior with  hyperparameters shape alpha_a and scale alpha_b.
#'  Precision parameters sigma_a,sigma_b, tau_a,tau_b are similar for the residual and penalties' precision, respectively.
#' 
#' @param formula Similar as for \code{\link[rsstap]{sstap_lmer}}, though fdp_staplmer is currently restricted to only one stap term.
#' @param benvo built environment - \code{\link[rbenvo]{benvo}} - object from containing the relevant subject - Built Environment data
#' @param weights weights for weighted regression - default is vector of ones 
#' @param alpha_a alpha gamma prior hyperparameter or alpha if fix_alpha = TRUE
#' @param alpha_b alpha gamma prior hyperparameter
#' @param sigma_a precision gamma prior hyperparameter
#' @param sigma_b precision gamma prior hyperparameter
#' @param tau_a penalty parameters gamma prior hyperparameter
#' @param tau_b penalty parameters gamma prior hyperparameter
#' @param K truncation number
#' @param iter_max maximum number of iterations
#' @param burn_in number of burn in iterations
#' @param thin number by which to thin samples
#' @param chains number of randomly initialized chains to run
#' @param fix_alpha boolean value indicating whether or not to fix the concentration parameter
#' @param seed random number generator seed will be set to default value if not by user
#' @param scale boolean determining if fixed effects matrix is scaled for estimation
#' @param center boolean determining if fixed effects matrix is centered for estimation
#' @param subsample_yhat  integer value indicating how many samples to subsample of yhat samples. Useful when N is big.
#' @param ... optional arguments to \code{\link{fdp_staplmer.fit}}
#' @importFrom stats is.empty.model model.matrix model.response as.formula gaussian terms
#' @export
#' @return a stapDP model object
#' 
fdp_staplmer <- function(formula,
						 benvo,
						 weights = NULL,
						 alpha_a = 1,
						 alpha_b = 1, 
						 sigma_a = 1,
						 sigma_b = 1,
						 tau_a = 1,
						 tau_b = 1,
						 K = 5L,
						 iter_max = 1000L,
						 burn_in = floor(iter_max/2),
						 thin = 1L,
						 chains = 1L,
						 fix_alpha = FALSE,
						 seed = NULL,
						 scale = TRUE,
						 center = TRUE,
						 subsample_yhat = NULL,
						 ...){

	## Parameter check
	num_posterior_samples <- sum((seq(from=burn_in+1,to=iter_max,by=1) %%thin)==0)
	stopifnot(burn_in<iter_max && burn_in >= 0)
	stopifnot(num_posterior_samples>0)
	stopifnot(all(c(alpha_a,alpha_b,sigma_a,sigma_b,tau_a,tau_b)>0))
	stopifnot(thin>0)
	## 
	call <- match.call(expand.dots = TRUE)
	
	spec <- get_stapDPspec(formula,K,benvo)
	foo <- spec$stapless_formula
	mf <- rbenvo::longitudinal_design(benvo,foo)
	W <- get_W(mf$glmod)
	subj_mat <- get_subjmat(mf$glmod)
	if(is.null(seed))
	  seed <- 3413431L
	if(is.null(subsample_yhat)){
		subsample_yhat <- 1:num_posterior_samples
		obs_ix <- 1:length(mf$y)
	}
	else {
	  if(length(subsample_yhat)>1){
	    obs_ix <- sample(1:length(mf$y),size = subsample_yhat[1],replace = F)
	    subsample_yhat <- subsample_yhat[2]
	  }
	  subsample_yhat <- sample(1:num_posterior_samples,size = subsample_yhat,replace=F)
	}
		

	if(vapply(list(mf$glmod$reTrms$flist),nlevels,1)>1)
		stop("Estimation of only 1 group term currently implimented")

	if(is.null(weights))
	  weights <- rep(1,length(mf$y))

	if(all(mf$X[,1]==1)){
		has_intercept <- TRUE
		Z <- scale(mf$X[,2:ncol(mf$X) , drop = F],scale = scale, center = center)
	}
	else{
		has_intercept <- FALSE 
		Z <- scale(mf$X, scale = scale, center = center)
	}
	Z_scl <- attr(Z,"scaled:scale")
	Z_cnt <- attr(Z,"scaled:center")
	if(is.null(Z_scl))
	  Z_scl <- rep(1,ncol(Z))
	if(is.null(Z_cnt))
	  Z_cnt <- rep(0,ncol(Z))
	if(has_intercept)
	  Z <- cbind(1,Z)

	rank <- spec$smooth_objs[[1]]$rank
	
	
	fit <- lapply(1L:chains,function(x) fdp_staplmer.fit(y = mf$y,
	                                                    Z = Z,
	                                                    X = spec$X,
	                                                    W = W,
	                                                    subj_mat = Matrix::t(subj_mat),
	                                                    subj_n = Matrix::rowSums(subj_mat),
	                                                    weights = weights, 
	                                                    alpha_a = alpha_a,
	                                                    alpha_b =  alpha_b,
	                                                    sigma_a = sigma_a,
	                                                    sigma_b = sigma_b,
	                                                    tau_a = tau_a,
	                                                    tau_b = tau_b,
	                                                    K = K,
														rank_one = rank[1],
														rank_two = rank[2],
	                                                    iter_max = iter_max,
	                                                    burn_in = burn_in,
	                                                    thin = thin,
	                                                    fix_alpha = fix_alpha,
	                                                    bw = has_bw(spec),
	                                                    seed = seed + x,
	                                                    chain = x,...))
	
    out <- lapply(fit,function(x) list(beta = x$beta,
							pi = x$pi,
							sigma = x$sigma,
							alpha = x$alpha,
						   scales_one = x$tau,
						   scales_two = x$tau2,
							tau_b  = if(has_bw(spec)) x$tau_b,
							tau_w = if(has_bw(spec)) x$tau_w,
							tau2_b = if(has_bw(spec)) x$tau2_b,
							tau2_w = if(has_bw(spec)) x$tau2_w,
							yhat = x$yhat[subsample_yhat,obs_ix,drop=F],
							subj_b = x$subj_b,
							subj_D = x$subj_D,
							pmat = x$PairwiseProbabilityMat,
							clabels = x$cluster_assignment
							))
	out <- list(call=call,
				pars=out,
				spec = spec,
				mf= mf,
				formula = formula,
				alpha_a = alpha_a,
				alpha_b = alpha_b,
				K = K,
				Z_scl = Z_scl,
				Z_cnt = Z_cnt,
				has_intercept = has_intercept
				)
	if(has_bw(spec)){
		lapply(1:chains,function(x) {
				   out$pars[[x]]$tau_b = fit[[x]]$tau_b
				   out$pars[[x]]$tau_w = fit[[x]]$tau_w
					 })
	}


	return(stapDP(out))
}


#'  Functional Dirichlet Process Spatial Temporal Aggregated Predictor Linear Mixed Effects Regression Model Fit
#' 
#' @param y vector of outcomes
#' @param Z design matrix
#' @param X stap design matrix
#' @param W group terms design matrix from \code{\link[lme4]{glFormula}}
#' @param subj_mat matrix indexing subject-measurement locations in (Z,X,W)
#' @param subj_n  vector of number of subject measurements
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param sigma_a precision gamma prior hyperparameter
#' @param sigma_b precision gamma prior hyperparameter
#' @param tau_a penalty parameters gamma prior hyperparameter
#' @param tau_b penalty parameters gamma prior hyperparameter
#' @param K truncation number for DP mixture components
#' @param rank_one rank of first smoothing matrix
#' @param rank_two rank of second smoothing matrix
#' @param threshold number of members per cluster at which cluster is included in regression
#' @param weights weights for weighted regression - default is vector of ones 
#' @param iter_max maximum number of iterations
#' @param burn_in number of iterations to burn-in
#' @param thin number by which to thin samples
#' @param fix_alpha boolean value 
#' @param logging boolean value indicating whether or not to do a sample iteration with diagnostic log messages
#' @param bw boolean value indicating whether or not subject decomposition is used
#' @param seed random number generator seed will be set to default value if not by user
#' @param chain chain label
#' @param logging boolean parameter indicating whether or not a single iteration should be run with print messages indicating successful completion of the Sampler's sub modules
#' @export
#' 
fdp_staplmer.fit <- function(y,Z,X,W,
                             subj_mat,
							 subj_n,
							 weights = rep(1,length(y)),
							 alpha_a = 1,
							 alpha_b = 1, 
							 sigma_a = 1,
							 sigma_b = 1,
							 tau_a = 1,
							 tau_b = 1,
							 K = 5L,
							 rank_one,
							 rank_two,
							 threshold = 0L,
							 iter_max,
							 burn_in,
							 thin = 1L,
							 fix_alpha = FALSE,
							 bw = FALSE,
							 seed = NULL,
							 chain = 1L,
							 logging = FALSE
							 ){

	stopifnot(c(sigma_a,sigma_b,tau_a,tau_b,alpha_a,alpha_b)>0)
	stopifnot(length(weights) == length(y))
  if(is.null(seed)){
    seed <- 3413
  }

  num_posterior_samples <- sum((seq(from=burn_in+1,to=iter_max,by=1) %%thin)==0)
  stopifnot(num_posterior_samples>0)

	if(bw){
	  X_b <- X[[1]]
	  X_w <- X[[2]]
	  fit <- stappDP_merdecomp(y = y ,Z = Z,
							   X_b = X_b, X_w = X_w,
							   W = W, w = weights,
							   subj_mat_ = subj_mat, subj_n = subj_n,
							   alpha_a =  alpha_a, alpha_b = alpha_b,
							   sigma_a = sigma_a,sigma_b = sigma_b,
							   tau_a = tau_a, tau_b = tau_b,
							   K = K, 
							   subset_one = rank_one, 
							   subset_two = rank_two,
							   threshold = threshold,
							   iter_max = iter_max, burn_in = burn_in,
							   thin = thin,seed = seed, 
							   num_posterior_samples = num_posterior_samples,
							   chain = chain, fix_alpha = fix_alpha,
							   logging = logging)

	}else{
	  fit <- stappDP_mer_fit(y = y,Z = Z, X = X, W = W,w = weights,
							 subj_mat_ = subj_mat,subj_n = subj_n,
							 alpha_a = alpha_a,alpha_b = alpha_b,
							 sigma_a = sigma_a, sigma_b = sigma_b,
							 tau_a = tau_a, tau_b = tau_b,
							 K = K, subset_one = rank_one, 
							 subset_two = rank_two,threshold = threshold,
							 iter_max = iter_max, burn_in = burn_in,
							 thin = thin,seed = seed, chain = chain,
							 num_posterior_samples = num_posterior_samples,
							 fix_alpha = fix_alpha,
							 logging = logging)
	}


  return(fit)
}
