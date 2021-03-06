% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ff_importutils.R
\name{get_sample}
\alias{get_sample}
\title{Read an FCS File and Do Some Processing}
\usage{
get_sample(
  fn,
  compensate = TRUE,
  transform = TRUE,
  derail = TRUE,
  nice.names = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{fn}{The fully qualified FCS filename}

\item{compensate}{Logical, apply the internal compensation matrix}

\item{transform}{Logical, apply linear transformation to scattering parameters
and biexponential transformation to fluorescence parameters}

\item{derail}{Logical, get rid of events on the FSC-A and SSC-A positive axes}

\item{nice.names}{Logical, swap 'desc' for 'name' in the flowFrame}

\item{verbose}{Logical, report progress on the console}
}
\value{
A flowFrame, properly processed.
}
\description{
Read a FCS file using \code{\link[flowCore]{read.FCS}} and return a
flowFrame.  Along the way, optionally apply compensation, transformation,
removal of 'rail' events, and apply some parameter renaming for convenience
and nice pictures.
}
\details{
This function is a convenience wrapper to apply commonly used processing
steps as FCS data are read in.  It's often important that all files in a project
get processed identically.  Consolidating most of the steps in one function helps
us avoid mistakes.
}
