% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_axes.R
\name{ax}
\alias{ax}
\title{Draw Pretty Axes (legacy base graphics)}
\usage{
ax(
  axis = 1,
  type = c("biexp", "asinh", "log", "linear", "none"),
  max_channel = 262143,
  ticksize = 2,
  ...
)
}
\arguments{
\item{axis}{Which axis to draw (1 = x, 2 = y)}

\item{type}{Type of transformation (one of \code{"biexp", "log", "linear"})}

\item{max_channel}{The maximum value that the flow cytometer can generate.  For most
BD instruments which have 18-bit A-D converters, that number is (2 ^ 8 - 1), or 262143}

\item{ticksize}{How big the tick marks should be (in relative coordinates)}

\item{...}{Additional graphical parameters to pass to \code{\link[graphics]{axis}}}
}
\description{
Draw axes for biexp, log or linear transformations.
}
\details{
This function is used by \code{\link{pplot}}, and only works with methods
in the standard R graphics package.  You can use it yourself
for figures you generate that use one of these transformations. To do so, first generate
a plot, and suppress the appropriate axis, e.g. \code{..., xaxt = 'n', ...}.
}
\examples{
# create an empty plot with suitable ranges
plot (0, 0, pch = '', xlim = c(0, 5.4), ylim = c(0, 5.4),
      xaxt = 'n', yaxt = 'n',
      xlab = 'biexponential', ylab = 'logarithmic')
# add the axes
ax(axis = 1, type = 'biexp')
ax(axis = 2, type = 'log')
}
