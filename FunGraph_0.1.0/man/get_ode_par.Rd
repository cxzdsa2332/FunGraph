% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODE_solving.R
\name{get_ode_par}
\alias{get_ode_par}
\title{calculation all LOP parameters in parallel mode}
\usage{
get_ode_par(data, times, order, reduction, parallel = TRUE)
}
\arguments{
\item{data}{matrix of observed data}

\item{times}{vector of time point}

\item{order}{scalar of LOP order}

\item{reduction}{use n/log(n) dimension reduction}

\item{parallel}{use parallel computation or not}
}
\value{
list contain variable selection results and LOP parameters for every row
}
\description{
calculation all LOP parameters in parallel mode
}
