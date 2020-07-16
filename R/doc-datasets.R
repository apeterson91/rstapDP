


#' Small dataset for use in \pkg{rstapDP} examples and vignettes.
#'
#' @name FFR_subjects
#'
#' @format  A data frame with 1000 rows and 3 columns
#' \describe{
#' \item{\code{id}}{:The subject unique identifier}
#' \item{\code{sex}}{:The measurement unique identifier}
#' \item{\code{BMI}}{:The Built Environment Unique identifier}
#' }
#' 
#' 
"FFR_subjects"


#' Small dataset for use in \pkg{rstapDP} examples and vignettes.
#'
#' @name FFR_distances
#' @format  A data frame with 9501 rows and 2 columns
#' \describe{
#' \item{ \code{id}}{: The subject unique identifier}
#' \item{\code{Distance}}{:The simulated distance between a hypothetical subject and fast food restaurant.}
#' }
#' 
#' 
"FFR_distances"

#' Samll benvo for use in \pkg{rstapDP} examples and vignettes.
#'
#' @name FFR_benvo 
#' @format  A benvo consisting of the \code{FFR_subjects} and \code{FFR_distances} dataframes
#' @details see \code{\link[rstapDP]{FFR_subjects}} and \code{\link[rstapDP]{FFR_distances}}
#' For information on the dataframes. See \code{\link[rbenvo]{Benvo}} for information on benvos.
#' 
#' 
"FFR_benvo"
