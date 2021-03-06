% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_mode.R
\name{assign_mode.stapDP}
\alias{assign_mode.stapDP}
\title{Estimates Posterior Mode Cluster Assignment}
\usage{
\method{assign_mode}{stapDP}(x, loss = "green", ...)
}
\arguments{
\item{x}{stapDP object}

\item{loss}{one of c("green","square") denoting which loss function should be used. Defaults to "green".}

\item{...}{optional arguments passed on to specific loss function}
}
\value{
vector of cluster assignments
}
\description{
Estimates Posterior Mode Cluster Assignment
}
