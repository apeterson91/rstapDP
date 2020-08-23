#' Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @details This function fits a linear model in a bayesian paradigm with
#' improper priors assigned to the "standard" regression covariates designated 
#' in the formula argument and a Dirichlet process prior with normal-gamma base measure 
#' assigned to the stap basis function expansion using penalized splines via \code{\link[mgcv]{jagam}}.
#'
#' The concentration parameter is assigned gamma prior with  hyperparameters shape alpha_a and scale alpha_b.
#'  Precision parameters sigma_a,sigma_b, tau_a,tau_b are similar for the residual and penalties' precision, respectively.
#' 
#' @param formula Similar as for \code{\link[rsstap]{sstap_lm}}, though fdp_staplm is currently restricted to only one stap term.
#' @param benvo built environment object from the rbenvo package containing the relevant data
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
#' @param fix_alpha boolean value indicating whether or not to fix the concentration parameter
#' @param seed random number generator seed will be set to default value if not by user
#' 
#' @importFrom stats is.empty.model model.matrix model.response as.formula gaussian terms
#' @export
#' @return a stapDP model object
#' 
fdp_staplm <- function(formula,
                       benvo,
                       weights = NULL,
					   alpha_a = 1,
					   alpha_b = 1, 
					   sigma_a = 1,
					   sigma_b = 1,
					   tau_a = 1,
					   tau_b = 1,
					   K = 5,
					   iter_max = 1E3,
					   burn_in = 5E2,
					   thin = 1,
					   fix_alpha = FALSE,
					   seed = NULL){

	## Parameter check
	stopifnot(burn_in<iter_max && burn_in > 0)
	stopifnot(any(c(alpha_a,alpha_b,sigma_a,sigma_b,tau_a,tau_b)>0))
	stopifnot(thin>0)
	## 
	
	spec <- get_stapDPspec(formula,K,benvo)
	foo <- spec$stapless_formula
	mf <- rbenvo::subject_design(benvo,foo)
	subj_mat <- get_subjmat(mf$glmod)

	if(length(spec$term)>1)
		stop("Only one stap/sap/tap term allowed")

	if(is.null(weights))
	  weights <- rep(1,length(mf$y))

	
	fit <- fdp_staplm.fit(y = mf$y,
						  Z = mf$X,
						  X=spec$X, 
						  S=spec$S,
						  weights,
						  alpha_a,
						  alpha_b,
						  sigma_a,
						  sigma_b,
	                      tau_a,tau_b,
						  K,iter_max,
						  burn_in,thin,
						  fix_alpha,seed)
	
    fit <- list(pars = list(beta = fit$beta,
							pi = fit$pi,
							sigma = fit$sigma,
							alpha = fit$alpha, 
							yhat = fit$yhat,
							cluster_mat = fit$cluster_assignment,
							scales =  fit$tau,
							pmat = fit$PairwiseProbabilityMat,
							clabels = fit$cluster_assignment
							),
				mf = mf,
				spec = spec,
				formula = formula,
				alpha_a = alpha_a,
				alpha_b = alpha_b,
				K = K
				)

	return(stapDP(fit))
}


#'  Functional Dirichlet Process Spatial Temporal Aggregated Predictor Linear Model Fit
#' 
#' @param y vector of outcomes
#' @param Z design matrix
#' @param X stap design matrix
#' @param S list of penalty matrices from \code{\link[mgcv]{jagam}} 
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
#' @param seed random number generator seed will be set to default value if not by user
#' @export
#' 
fdp_staplm.fit <- function(y,Z,X,S,
                           weights = rep(1,length(y)),
                           alpha_a = 1,
                           alpha_b = 1, 
                           sigma_a = 1,
                           sigma_b = 1,
                           tau_a = 1,
                           tau_b = 1,
                           K = 5,
                           iter_max,
						   burn_in,
                           thin = 1,
						   fix_alpha = FALSE,
						   seed = NULL){

	stopifnot(c(sigma_a,sigma_b,tau_a,tau_b,alpha_a,alpha_b)>0)
	stopifnot(nrow(S[[1]]) == ncol(Z) + ncol(X[[1]])*K)
	stopifnot(length(weights) == length(y))
  if(is.null(seed)){
    seed <- 3413
  }

  num_posterior_samples <- sum((seq(from=burn_in+1,to=iter_max,by=1) %%thin)==0)
  stopifnot(num_posterior_samples>0)

  num_penalties <- length(S) ## default for smoothing
  S <- do.call(cbind,S)
  X <- do.call(cbind,X)
  fit <- stappDP_fit(y,Z,X,S,weights,alpha_a,alpha_b,
					 sigma_a,sigma_b,tau_a,tau_b,
					 K,num_penalties,iter_max,burn_in,
					 thin,seed,num_posterior_samples,fix_alpha)


  return(fit)
}
