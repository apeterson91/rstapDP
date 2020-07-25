#' get_stapless_formula
#'
#' Get formula for typical covariates
#'
#' @keywords internal
#' @param f formula from stap_glm
#' @return formula without  stap(),sap(),tap() components
#'
get_stapless_formula <- function(f){
    
    with_bars <- f
    f <- lme4::nobars(f)
	get_ics <- function(f,vec_var){
		which(all.names(f) %in% vec_var)
	}
    stap_ics <- get_ics(f, c("stap","stap_bw"))
    sap_ics <- get_ics(f,c("sap","sap_bw"))
    tap_ics <- get_ics(f,c("tap","tap_bw"))
    if(!length(stap_ics) & !length(sap_ics) & !length(tap_ics))
        stop("No covariates designated as 'stap','sap',or 'tap'  in formula", .call = F)
    stap_nms <- all.names(f)[stap_ics + 1]
	stap_bw <- (all.names(f)[stap_ics] %in% c("stap_bw"))*1
    sap_nms <- all.names(f)[sap_ics + 1]
	sap_bw <- (all.names(f)[sap_ics] %in% c("sap_bw"))*1
    tap_nms <- all.names(f)[tap_ics + 1]
	tap_bw <- (all.names(f)[tap_ics] %in% c("tap_bw"))*1
	if(length(stap_nms)>0){
		stap_nms <- cbind(stap_nms,"Distance-Time",stap_bw)
	}
	if(length(sap_nms)>0)
		sap_nms <- cbind(sap_nms,"Distance",sap_bw)
	if(length(tap_nms)>0)
		tap_nms <- cbind(tap_nms,"Time",tap_bw)

	stap_mat <-rbind(stap_nms,sap_nms,tap_nms)

    not_needed <- c(stap_nms,sap_nms,tap_nms)
    formula_components <- all.vars(f)[!(all.vars(f) %in% not_needed)]
    if(!attr(terms(f),"intercept"))
        formula_components <- c(formula_components,"0")
    if(grepl("cbind",all.names(f))[2]){
        new_f1 <- paste0("cbind(",formula_components[1],", ",formula_components[2], ")", " ~ ")
        ix <- 3
    }
    else{
        new_f1 <- paste0(formula_components[1],' ~ ')
        ix <- 2
    }

    new_f2 <- paste(formula_components[ix:length(formula_components)],collapse = "+")
    new_f <- paste0(new_f1,new_f2)

	str <- purrr::map(stap_mat[,2],function(x) {
	  switch(x,
	         "Distance-Time"= "t2(Distance,Time,bs='ps')", 
	         "Distance" = "s(Distance,bs='ps',k=5)",
	         "Time"= "s(Time,bs='ps')")
	  })

	fake_formula <- purrr::map(str,function(x) as.formula(paste0("ID~ -1 + ",paste0(x,collapse="+"))))

    return(
		   list(stapless_formula = as.formula(new_f, env = environment(f)),
				fake_formula = fake_formula,
				stap_mat = stap_mat
			   )
		   )
}
