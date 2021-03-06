% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_histogram_analysis.R
\name{find.local.minima}
\alias{find.local.minima}
\title{Find Local Minima in a KDE}
\usage{
find.local.minima(kde, thresh = 0.001, show = FALSE, ...)
}
\arguments{
\item{kde}{The input Kernel Density Estimate}

\item{thresh}{Threshold above which valleys are considered}

\item{show}{Logical, should we display our work?}

\item{...}{Additional arguments passed to \code{par() \link[graphics]{par}}, Relevant only if show = TRUE.}
}
\value{
A list with three elements:
\describe{
\item{x}{X coordinates of the detected peaks}
\item{y}{Y coordinates of the detected peaks}
\item{pick}{The indices in the original kde corresponding to the peaks}
}
}
\description{
Find the valleys in a Kernel Density Estimate.
}
