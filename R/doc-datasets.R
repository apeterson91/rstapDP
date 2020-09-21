


#' Small dataset for use in \pkg{rstapDP} examples and vignettes.
#'
#' @name FFR_subjects
#'
#' @format  A data frame with 1000 rows and 3 columns
#' \describe{
#' \item{\code{id}}{The subject unique identifier}
#' \item{\code{sex}}{The measurement unique identifier}
#' \item{\code{BMI}}{The Built Environment Unique identifier}
#' }
#' 
#' 
"FFR_subjects"


#' Small dataset for use in \pkg{rstapDP} examples and vignettes.
#'
#' @name FFR_distances
#' @format  A data frame with 9501 rows and 2 columns
#' \describe{
#' \item{ \code{id}}{ The subject unique identifier}
#' \item{\code{Distance}}{The simulated distance between a hypothetical subject and fast food restaurant.}
#' }
#' 
#' 
"FFR_distances"

#' Samll benvo for use in \pkg{rstapDP} examples and vignettes.
#'
#' @name FFR_benvo 
#' @format  A \code{\link[rbenvo]{benvo}} consisting of the \code{FFR_subjects} and \code{FFR_distances} dataframes
#' @details see \code{\link{FFR_subjects}} and \code{\link{FFR_distances}}
#' For information on the dataframes.
#' @seealso  \code{\link[rbenvo]{benvo}} for information on benvos.
#' 
#' 
"FFR_benvo"

#' Longitudinal Clustering benvo for use in \pkg{rstapDP} vignettes.
#'
#' @name longitudinal_clusters
#' @format A \code{\link[rbenvo]{benvo}} consisting of a subject and subject-distance dataframes
#' @details The outcome in the subject dataframe is simulated using a cluster specific effect as follows:
#' \deqn{E[Y_{ij}|b_{i1},b_{i2}] = 33  - 2.2 \text{sex}_i + .1*(year_{ij}) +  f_k(\text{FFR exposure}) +b_{i1} + b_{i2}(year_{ij})  } Where \eqn{k= 1,2,3} with varying probabilities.
#' @seealso the generating \href{https://github.com/apeterson91/rstapDP/tree/master/data-raw}{code} on Github.
#'
"longitudinal_clusters"


#' Complex Longitudinal Clustering benvo for use in \pkg{rstapDP} vignettes.
#'
#' @name complex_longitudinal_clusters
#' @format A \code{\link[rbenvo]{benvo}} consisting of a subject and subject-distance dataframes
#' @details The outcome in the subject dataframe is simulated using a cluster specific effect as follows:
#' \deqn{E[Y_{ij}|b_{i1},b_{i2}] = 33  - 2.2 \text{sex}_i + .1*(year_{ij}) +  \bar{f}_k(FFR_i) +  \Delta f_k(FFR_{ij}) +b_{i1} + b_{i2}(year_{ij})  } Where \eqn{k= 1,2,3} with varying probabilities.
#' @seealso the generating \href{https://github.com/apeterson91/rstapDP/tree/master/data-raw}{code} on Github, as well as the Longitudinal vignette.
#'
"complex_longitudinal_clusters"
