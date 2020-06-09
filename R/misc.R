#' smooth function for stap terms
#' 
#' 
#' 
# s_stap <- function(x,knots){
#   return(list(term = mgcv::s(x)$term))
#   # Not using right now
#   #mgcv::s(x,k=-1,fx=FALSE,bs=bs,m=NA,by=NA,xt=NULL,id=NULL,sp=NULL,pc=NULL)
#}

#' Retrieves the design matrix
#'
#'
get_X <- function(dt_data,subject_data,knots){

  dt <- subject_data %>% dplyr::left_join(dt_data)
  
  X <- dt %>% dplyr::group_by(id) %>% 
    dplyr::summarise(D_int = sum(!is.na(Distance)),
              Linear = sum(Distance,na.rm=T)) %>% dplyr::select(D_int,Linear) %>% 
    as.matrix()
  if(all(X[,1]==1)) ## identifiability constraint
    X <- X[,2,drop=F]
  
  X_<- dt %>% tidyr::crossing(dplyr::tibble(knot = knots,
                                            k_id=1:length(knots))) %>% 
    dplyr::group_by(id,k_id) %>% 
    dplyr::mutate(Intercept = n(),
                  Linear = sum(Distance)) %>% 
    dplyr::summarise(Exposure = sum(abs(Distance - knot)^3)) %>% 
    dplyr::ungroup() %>% 
    tidyr::spread(k_id,Exposure) %>% dplyr::select(-id) %>% 
    as.matrix()
  X <- cbind(X,X_)
  
  stopifnot(nrow(X)==nrow(subject_data))
  return(X)
	# if(is.null(dt_data))
	# 	X <- smooth.construct(stap_term,subject_data,knots=NULL)$X
	
}

get_predict_X <- function(d,knots){
  mat <- cbind(1,d,outer(d,knots,FUN = function(x,y) abs(x-y)^3))
}

get_D <- function(knots){
  D <- outer(knts,knts,FUN = function(x,y) abs(x-y)^3) + diag(length(knots))
}

get_X_transformed <- function(X,D){
  p <- ncol(D)
  pp <- ncol(X)
  start <- pp-p+1
  X_t <- cbind(X[,1:(start-1)],(X[,start:pp] %*% t(chol(solve(D))) ))
}

transform_beta <- function(beta,D){
  start <- ncol(beta) - ncol(D) + 1
	bt <- cbind(beta[,1:(start-1)],t(t(chol(D)) %*% t(beta[,start:ncol(beta)])))
}

