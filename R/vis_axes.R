# vis_axes.R
#
# generates nicely labeled axes.
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#

#' @title Draw Pretty Axes
#' @description Draw axes for biexp, log or linear transformations.
#' @param axis Which axis to draw (1 = x, 2 = y)
#' @param type Type of transformation (one of \code{"biexp", "log", "linear"})
#' @param max_channel The maximum value that the flow cytometer can generate.  For most
#' BD instruments which have 18-bit A-D converters, that number is (2 ^ 8 - 1), or 262143
#' @param ticksize How big the tick marks should be (in relative coordinates)
#' @param ... Additional graphical parameters to pass to \code{\link[graphics]{axis}}
#' @details This function is used by \code{\link{pplot}}, and only works with methods
#' in the standard R graphics package.  You can use it yourself
#' for figures you generate that use one of these transformations. To do so, first generate
#' a plot, and suppress the appropriate axis, e.g. \code{..., xaxt = 'n', ...}.
#' @examples
#' # create an empty plot with suitable ranges
#' plot (0, 0, pch = '', xlim = c(0, 5.4), ylim = c(0, 5.4),
#'       xaxt = 'n', yaxt = 'n',
#'       xlab = 'biexponential', ylab = 'logarithmic')
#' # add the axes
#' ax(axis = 1, type = 'biexp')
#' ax(axis = 2, type = 'log')
#' @export
ax <- function (axis = 1, type = c("biexp", "asinh", "log", "linear"), max_channel = 262143, ticksize = 2, ...) {

  type       = match.arg(type)

  channels = max_channel
  decades=8

  all.ticks<-NULL

  if (!(type %in% c("biexp", "asinh", "log", "linear"))) {
    cat ("invalid axis type\n")
    return()
  }
  if (type == "log") {
    start_decade=0
    major<-(10^(start_decade:decades))
    for(i in 1:(length(major)-1)){
      all.ticks<-c(all.ticks,seq(major[i],major[i+1],l=10))
    }
    all.ticks <- log(all.ticks, base=10); #Log base 10
    major     <- log(major, base=10);
  }
  if (type == "biexp") {
    start_decade=2
    major<-(10^(start_decade:decades))
    for(i in 1:(length(major)-1)){
      all.ticks<-c(all.ticks,seq(major[i],major[i+1],l=10))
    }
    all.ticks <- biexp.transform(all.ticks);
    major     <- biexp.transform(major);
  }
  if (type == "asinh") {
    start_decade=2
    major<-(10^(start_decade:decades))
    for(i in 1:(length(major)-1)){
      all.ticks<-c(all.ticks,seq(major[i],major[i+1],l=10))
    }
    all.ticks <- asinh.transform(all.ticks);
    major     <- asinh.transform(major);
  }
  if (type == "linear") {
    axis (side=axis)
    return()
  }
  axis(side=axis, all.ticks, label=FALSE, tcl=-0.25)
  axis(side=axis, at=major, lwd.ticks=ticksize, label=tick.labels(start_decade, decades),...)
}


#' @title Label (Quasi-)Logarithmic Tickmarks
#' @description Create powers-of-ten tickmark labels
#' @param first_decade First decade (e.g. 1 would correspond to 10^1)
#' @param last_decade Last decade (e.g. 5 would correspond to 10^5)
#' @return A list of parsed labels suitable for the "label" parameter of \link{ax}
#' or \link[graphics]{axis}.
#' @export
tick.labels <- function (first_decade, last_decade) {
  label <- vector('character')
  n_labels <- last_decade - first_decade + 1
  for (i in 1:n_labels) {
    label[i] <- paste ("10^", first_decade + i - 1, sep="")
  }
  ex = parse(text=label)

  return(ex)
}
