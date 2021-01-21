#' Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Mixed Effects Regression Model
#' 
#' @details This function fits a linear mixed effects regression model in a Bayesian paradigm with
#' improper priors assigned to the "standard" regression covariates designated 
#' in the formula argument and a Dirichlet process prior with normal-gamma base measure 
#' assigned to the stap basis function expansion using penalized splines via \code{\link[mgcv]{jagam}}.
#' normal priors are placed on the latent group variables and an improper prior is placed on the 
#' correlation matrix leading to a Wishart posterior.
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
						 ...){

	## Parameter check
	stopifnot(burn_in<iter_max && burn_in > 0)
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


	if(vapply(list(mf$glmod$reTrms$flist),nlevels,1)>1)
		stop("Estimation of only 1 group term currently implimented")

	if(is.null(weights))
	  weights <- rep(1,length(mf$y))

	if(all(mf$X[,1]==1)){
		has_intercept <- TRUE
		Z <- cbind(1,scale(mf$X[,2:ncol(mf$X) , drop = F],scale = scale, center = center))
	}
	else{
		has_intercept <- FALSE 
		Z <- scale(mf$X, scale = scale, center = center)
	}
	Z_scl <- attr(Z,"scaled:scale")
	Z_cnt <- attr(Z,"scaled:center")
	Z_scl <- ifelse(is.null(Z_scl),1,Z_scl)
	Z_cnt <- ifelse(is.null(Z_cnt),0,Z_cnt)
	
	
	fit <- lapply(1L:chains,function(x) fdp_staplmer.fit(y = mf$y,
	                                                    Z = Z,
	                                                    X = spec$X,
	                                                    W = W,
	                                                    S = spec$S,
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
							scales = x$tau,
							tau_b  = if(has_bw(spec)) x$tau_b,
							tau_w = if(has_bw(spec)) x$tau_w,
							yhat = x$yhat,
							subj_b = x$subj_b,
							subj_D = x$subj_D,
							pmat = x$PairwiseProbabilityMat,
							clabels = x$cluster_assignment
							))
	out <- list(pars=out,
				spec = spec,
				mf= mf,
				formula = formula,
				alpha_a = alpha_a,
				alpha_b = alpha_b,
				benvo = benvo,
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
#' @param S list of penalty matrices from \code{\link[mgcv]{jagam}} 
#' @param subj_mat matrix indexing subject-measurement locations in (Z,X,W)
#' @param subj_n  vector of number of subject measurements
#' @param alpha_a alpha gamma prior hyperparameter
#' @param alpha_b alpha gamma prior hyperparameter
#' @param sigma_a precision gamma prior hyperparameter
#' @param sigma_b precision gamma prior hyperparameter
#' @param tau_a penalty parameters gamma prior hyperparameter
#' @param tau_b penalty parameters gamma prior hyperparameter
#' @param K truncation number for DP mixture components
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
#' @param summarize_yhat boolean value indicating whether a single mean vector of yhat values should be returned instead of a N X num samples matrix. Useful in situations where N is large.
#' @export
#' 
fdp_staplmer.fit <- function(y,Z,X,W,S,
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
              							 iter_max,
              							 burn_in,
              							 thin = 1L,
              							 fix_alpha = FALSE,
              							 bw = FALSE,
              							 seed = NULL,
              							 chain = 1L,
              							 logging = FALSE,
              							 summarize_yhat = FALSE
              							 ){

	stopifnot(c(sigma_a,sigma_b,tau_a,tau_b,alpha_a,alpha_b)>0)
	stopifnot(length(weights) == length(y))
  if(is.null(seed)){
    seed <- 3413
  }

  num_posterior_samples <- sum((seq(from=burn_in+1,to=iter_max,by=1) %%thin)==0)
  stopifnot(num_posterior_samples>0)

	if(bw){
	  num_penalties <- length(S[[1]])
	  S_b <- do.call(cbind,S[[1]])
	  S_w <- do.call(cbind,S[[2]])
	  X_b <- X[[1]]
	  X_w <- X[[2]]
	  fit <- stappDP_merdecomp(y,Z,X_b,X_w,W,S_b,S_w,
							   weights,subj_mat,
							   subj_n,alpha_a,alpha_b,
							   sigma_a,sigma_b,tau_a,tau_b,
							   K,num_penalties,iter_max,burn_in,
							   thin,seed,num_posterior_samples,chain,fix_alpha
							  )

	}else{
	  num_penalties <- length(S) ## default for smoothing
	  S <- do.call(cbind,S)
	  X <- do.call(cbind,X)
	  fit <- stappDP_mer_fit(y,Z,X,W,S,weights,subj_mat,
							 subj_n,alpha_a,alpha_b,
							 sigma_a,sigma_b,tau_a,tau_b,
							 K,num_penalties,iter_max,burn_in,
							 thin,seed,chain,num_posterior_samples,
							 fix_alpha,logging,summarize_yhat)
	}


  return(fit)
}
