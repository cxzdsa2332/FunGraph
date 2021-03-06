% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODE_solving.R
\name{LOP_rk4}
\alias{LOP_rk4}
\title{generate LOP based rk4 with step=h}
\usage{
LOP_rk4(x0, y0, h, d_LOP_fit)
}
\arguments{
\item{x0}{scalar of x (equal tp a time point)}

\item{y0}{scalar of an observation}

\item{h}{scalar of step h(x1-x0)}

\item{d_LOP_fit}{derivative of LOP}
}
\value{
scalar of y value (y next)
}
\description{
generate LOP based rk4 with step=h
}
