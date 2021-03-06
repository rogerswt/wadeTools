% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans_asinh.R
\name{w.arcsinh}
\alias{w.arcsinh}
\title{Calculate a Arcsinh Transformation}
\usage{
w.arcsinh(ff, cofactor = 5, params)
}
\arguments{
\item{ff}{a flowFrame}

\item{params}{The parameters of ff to transform}

\item{a}{The "a" parameter.  Do not override unless you know what you're doing.}
}
\value{
A transformed flowFrame
}
\description{
Do the arcsinh mathematical transform

This function performs the arcsinh transform directly on ff.
See \code{\link{biexpTransform}} for flowCore compatibility.
}
