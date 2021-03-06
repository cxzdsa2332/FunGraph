% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biFunClu.R
\name{get_cluster}
\alias{get_cluster}
\title{biFunClu function}
\usage{
get_cluster(input, itermax = 100, Delta = 100, write_data = TRUE)
}
\arguments{
\item{input}{list object from get_init_par function contain initial parameters and data}

\item{itermax}{scalar constrain the maximum number of iteracion}

\item{Delta}{scalar prevent Log-likelihood drastically change}

\item{write_data}{write result table for further analysis}
}
\value{
Log-likelihood value
}
\description{
biFunClu function
}
