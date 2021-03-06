% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biFunClu.R
\name{get_BIC}
\alias{get_BIC}
\title{batch biFunClu}
\usage{
get_BIC(data, order, times, rep, min_cluster, max_cluster, itermax)
}
\arguments{
\item{data}{dataframe contain numeric values for further analysis}

\item{order}{scalar of LOP order}

\item{times}{vector of time points}

\item{rep}{numeric value for repetation number under a k}

\item{min_cluster}{scalar of minimum cluster to calculate(must >= 2)}

\item{max_cluster}{scalar of maximum cluster to calculate}

\item{itermax}{scalar constrain the maximum number of iteracion}
}
\value{
batch biFunClu results in a list
}
\description{
batch biFunClu
}
