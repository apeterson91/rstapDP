
#' Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @details This function fits a linear model in a bayesian paradigm with
#' improper priors assigned to the "standard" regression covariates designated 
#' in the formula argument and a Dirichlet process prior with normal base measure and 
#' scale tau_df if penalize = F or inverse-chi square(tau_df) otherwise assigned to the 
#' stap_term basis function expansion. 
#' The concentration parameter is assigned the hyperparameters alpha_a and alpha_b.
#' The residual variance is also assigned an improper prior.
#' 
#' @param formula Similar as for \code{\link[rsstap]{sstap_lm}}. 
#' @param benvo built environment object from the rbenvo package containing the relevant data
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
#' @return a stapDP model object
#' 
fdp_staplm <- function(formula,
					   benvo,
					   weights = NULL,
					   tau_df = 1,
					   alpha_a = 1,
					   alpha_b = 1, 
					   K = 5,
					   penalize = T,
					   iter_max = 1E3,
					   burn_in = 5E2,
					   thin = 1,
					   id_colname = 'id',
					   seed = NULL){

	## Parameter check
	stopifnot(burn_in<iter_max && burn_in > 0)
	stopifnot(alpha_a>0)
	stopifnot(alpha_b>0)
	stopifnot(thin>0)
	## 
	
	foo <- get_stapless_formula(formula)
	f <- foo$stapless_formula
	mf <- rbenvo::subject_design(benvo,f)
	Z <- mf$X
	call <- match.call(expand.dots = TRUE)
	if(nrow(foo$stap_mat)>1)
		stop("Only one stap/sap/tap term allowed")
	stap_term <- foo$stap_mat[,1]
	stap_component <- foo$stap_mat[,2]
	bw <- as.integer(foo$stap_mat[,3])
	stap_formula <- foo$fake_formula[[1]]

	
	## Handle Zero exposure in dt_data
	jd <- mgcv::jagam(formula = stap_formula,family = gaussian(),
					  data = rbenvo::joinvo(benvo,
											stap_term,
											stap_component,
											NA_to_zero = TRUE),
					  file = tempfile(fileext=".jags"),
					  weights = NULL,
					  offset = NULL,
					  centred = FALSE,
					  diagonalize = FALSE)

	X <- rbenvo::aggrenvo(benvo,jd$jags.data$X,stap_term,stap_component)
	S <- lapply(jd$pregam$S,function(x) kronecker(diag(K),x))
	S <- lapply(S,function(m) rbind(matrix(0,ncol=ncol(m)+ncol(Z),nrow=ncol(Z)),
	           cbind(matrix(0,nrow=nrow(m),ncol=ncol(Z)),m)))
	
	if(is.null(weights))
	  weights <- rep(1,length(mf$y))
	
	
	fit <- fdp_staplm.fit(y = mf$y,Z,X, S, alpha_a,alpha_b,K,penalize,tau_df,weights,iter_max,burn_in,thin,seed)
	
    fit <- list(beta = fit$beta,
                probs = fit$pi,
                sigma = fit$sigma,
                alpha = fit$alpha,
                yhat = fit$yhat,
                cluster_mat = fit$cluster_assignment,
                scales = if(penalize) fit$tau else matrix(rep(tau_df,K),ncol=K),
				num_penalties = length(S),
                pmat = fit$PairwiseProbabilityMat,
                clabels = fit$cluster_assignment,
                Znames = colnames(Z),
                ncol_X = ncol(X),
                penalize = penalize,
                formula = formula,
                sobj = jd$pregam$sobj,
        				alpha_a = alpha_a,
        				alpha_b = alpha_b,
                y = mf$y,
                K = K)
	return(stapDP(fit))
}


#'  Functional Dirichlet Process Spatial Temporal Aggregated Predictor in a Linear Model
#' 
#' @param y vector of outcomes
#' @param Z design matrix
#' @param X stap design matrix
#' @param S list of penalty matrices from \code{\link[mgcv]{jagam}} 
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
	stopifnot(nrow(S) == ncol(Z) + ncol(X)*K)
	stopifnot(length(w) == length(y))
  if(is.null(seed)){
    seed <- 3413
  }

  num_posterior_samples <- sum((seq(from=burn_in+1,to=iter_max,by=1) %%thin)==0)

  if(penalize){
	  num_penalties <- length(S) ## default for smoothing
	  S <- do.call(cbind,S)
	  fit <- stappDP_fit(y,Z,X,S,w,tau_df,alpha_a,alpha_b,K,num_penalties,iter_max,burn_in,thin,seed,num_posterior_samples)
  }else
	  fit <- stapDP_fit(y,Z,X,w,
						  tau_df,alpha_a,alpha_b,
						  K,iter_max,burn_in,thin,
						  seed,num_posterior_samples)


  return(fit)
}
