# histogram_analysis.R
#
# Utility functions for dealing with histograms in the form of
# Kernel Density Estimates (kde's).
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#


#' @title Find Local Maxima in a KDE
#' @description Find the peaks in a Kernel Density Estimate.
#' @param kde The input Kernel Density Estimate
#' @param thresh Threshold above which peaks are considered
#' @param show Logical, should we display our work?
#' #' @param ... Additional arguments passed to par().  Relevant only if show = TRUE.
#' @return A list with three elements:
#' \describe{
#'   \item{x}{X coordinates of the detected peaks}
#'   \item{y}{Y coordinates of the detected peaks}
#'   \item{pick}{The indices in the original kde corresponding to the peaks}
#' }
#' @export
find.local.maxima <- function (kde, thresh=.05, show=FALSE, ...) {

  ind <- msExtrema(kde$y, span=11)$index.max

  max_y <- max (kde$y)
  sel <- kde$y / max_y > thresh

  ind <- which (ind & sel)

  if (show) {
    par(...)
    plot (kde, type='l', xlab="", ylab="", xaxt='n')
    points (kde$x[ind], kde$y[ind], pch=20, cex=.5, col='red')
  }

  return (list(x=kde$x[ind], y=kde$y[ind], pick=ind))
}

#' @title Find Local Minima in a KDE
#' @description Find the valleys in a Kernel Density Estimate.
#' @param kde The input Kernel Density Estimate
#' @param thresh Threshold above which valleys are considered
#' @param show Logical, should we display our work?
#' @param ... Additional arguments passed to par().  Relevant only if show = TRUE.
#' @return A list with three elements:
#' \describe{
#'   \item{x}{X coordinates of the detected peaks}
#'   \item{y}{Y coordinates of the detected peaks}
#'   \item{pick}{The indices in the original kde corresponding to the peaks}
#' }
#' @export
find.local.minima <- function (kde, thresh=.001, show=FALSE, ...) {

  ind <- msExtrema(kde$y, span=11)$index.min

  max_y <- max (kde$y)
  sel <- kde$y / max_y > thresh

  pick <- which (ind & sel)

  if (show) {
    par(...)
    plot (kde, type='l', xlab="", ylab="", xaxt='n')
    points (kde$x[pick], kde$y[pick], pch=20, cex=.5, col='red')
  }

  return (list(x=kde$x[pick], y=kde$y[pick], pick=pick))
}

# helper function, not to be exposed
msExtrema <- function(x, span=3) {
  # WTR 2016-08-26  - looks like peaks has been deprecated again!
  # Found a function in a previous version and copied it here.  It depends on
  # splus2r
  # find local maxima
  index1 <- peaks(x, span=span, strict=FALSE)

  # find local minima
  index2 <- peaks(-x, span=span, strict=FALSE)

  # remove the interior of plateaus
  index.max <- index1 & !index2
  index.min <- index2 & !index1

  # construct output
  list(index.max=index.max, index.min=index.min)
}


#' @title Calculate the Width of a KDE
#' @description  find the full width half-maximum (or other width at percent maximum) of a kde.
#' @param kde The input kde
#' @param frac The height ([0, 1]) at which to compute the width.
#' @return A list with 4 elements:
#' \describe{
#'   \item{center}{X coordinate of the peak center}
#'   \item{left}{X coordinate of the left intercept}
#'   \item{right}{X coordinate of the right intercept}
#'   \item{frac}{The value of the input parameter 'frac'}
#' }
#' @export
find.width = function (kde, frac=0.5) {
  x = kde$x
  y = kde$y
  mx = max(kde$y)

  for (i in 1:length(y)) {
    if (y[i]/mx > frac) {break}
  }
  left = x[i]
  for (i in length(y):1) {
    if (y[i]/mx > frac) {break}
  }
  right = x[i]
  center = left + (right - left)/2

  return (list(center=center, left=left, right=right, frac=frac))
}

#' @title Calculate the Cummulative Probability Distribution of a KDE
#' @description The cummulative distribution function (think of the kde as a
#' differential probability distribuion).
#' @param kde The input kde
#' @return A list with 2 elements:
#' \describe{
#'   \item{x}{X coordinates of the CDF}
#'   \item{y}{Y coordinates of the CDF, normalized to 1.0.}
#' }
#' @export
cumm.kde = function (kde) {
  n = length(kde$y)
  cval = vector('numeric', length=n)
  cval[1] = kde$y[1]
  for (i in 2:n) {
    cval[i] = cval[i-1] + kde$y[i]
  }
  cval = cval / cval[n]

  return(list(x=kde$x, y=cval))
}

#' @title Estimate the First Derivative of a KDE
#' @description Estimate the first derivative of a kde using the symmetric difference quotient.
#' @param kde The input kde
#' @param normalize Logical, should the resulting kde be normalized to 1.0?
#' @return Another kde.
#' @export
deriv.kde = function (kde, normalize=TRUE) {
  x = kde$x
  y = kde$y
  npts = length(y)
  yp = vector('numeric')
  for (i in 2:(npts-1)) {
    yp[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
  }
  yp[1] = yp[2]
  yp[npts] = yp[npts-1]

  if (normalize) {
    yp = yp / max (abs(yp))
  }
  res = list(x=x, y=yp)
  res
}

#' @title Normalize a KDE
#' @description Normalize a kde to 1.0
#' @param kde The input kde
#' @return A normalized version of the kde.
#' @export
normalize.kde = function(kde) {
  kde$y = kde$y / max(kde$y)

  kde
}
