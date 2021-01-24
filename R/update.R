
#' Update STAPDP Model fit
#' 
#' To save time from creating the data structures that are used in stapDP
#' models, one can "re-use" the previous data structure but change aspects of the 
#' fixed regression coefficients, hyperparameter or sampler settings.
#' 
#' @export
#' @param object stapDP object
#' @param formula formula used to specific fixed effects, see e.g. \code{\link[stats]{lm}}
#' @param newdata subject dataframe for subject level effects
#'
update.stapDP <- function(object,
						  formula , 
						  newdata,
						  weights = NULL,
						  iter_max = 1E3,
						  burn_in = floor(iter_max/2),
						  chains = 1,
						  alpha_a = 1,
						  alpha_b = 1,
						  sigma_a = 1,
						  sigma_b = 1,
						  tau_a = 1,
						  tau_b = 1,
						  K = 5,
						  thin = 1,
						  fix_alpha = FALSE,
						  seed = NULL,
						  scale = TRUE,
						  center = TRUE,
						  ... ){

	stopifnot(burn_in<iter_max && burn_in > 0)
	stopifnot(all(c(alpha_a,alpha_b,sigma_a,sigma_b,tau_a,tau_b)>0))
	stopifnot(thin>0)
	spec <- object$spec

	if(is.null(newdata))
		newdata <- object$spec$mf$X
	if(is.null(seed))
		seed <- 3413431L
	if(is.null(weights))
		weights <- rep(1,nrow(newdata))

	X <- object$spec$X
	if(is.mer(object)){
		mf <- .design(formula,newdata)
		spec$mf <- mf
		W <- get_W(mf$glmod)
		subj_mat <- get_subjmat(mf$glmod)
		Z <- mf$X
		y <- mf$y
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

		fit <- lapply(1L:chains,function(x) fdp_staplmer.fit(y = y,
															Z = Z,
															X = X,
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
															iter_max = iter_max,
															burn_in = burn_in,
															thin = thin,
															fix_alpha = fix_alpha,
															bw = has_bw(spec),
															seed = seed + x,
															chain = x,...))
	}
	else{
		mf <- model.frame(formula,subj_newdata)
		##TODO: Impliment update for single measurement newdata
	}

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
				K = K,
				Z_scl = Z_scl,
				Z_cnt = Z_cnt,
				has_intercept = has_intercept
				)
	return(stapDP(out))

}

.design <- function(formula,newdata){

	  control <- glmerControl(check.nlev.gtreq.5 = "ignore",
	                             check.nlev.gtr.1 = "stop",
	                             check.nobs.vs.rankZ = "ignore",
	                             check.nobs.vs.nlev = "ignore",
	                             check.nobs.vs.nRE = "ignore" )

	  mf <- lme4::glFormula(formula = formula,
							data = newdata, control = control)
	  y <- mf$fr[,as.character(mf$formula[2L])]
	  X <-  mf$X
	  out <- list(y=y,X=X,glmod=mf)
	}
