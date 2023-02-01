# vis_pplot.R
#
#	This function makes a "FlowJo" style plot.
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without Still Pond Cytomics express written consent.             ##
################################################################################
################################################################################
#
#	This function is a wrapper function for flowCore:plot
#	It creates a plot that looks very much like FlowJo.
#
#' @title Draw a FlowJo-style (sort of) Plot
#' @description Draw a picture, with colors reminiscent of FlowJo.  \strong{NOTE: this
#' function is deprecated.  Please consider using \link[wadeTools]{ggflow} instead.}
#' @param ff The flowFrame to be plotted
#' @param parameters A vector of two parameters in ff.  May be either character or numeric
#' @param blueBackground Logical, should the background be blue.  If false, you'll get a white
#' background (like FlowJo)
#' @param showZero Logical, draw dot-dash lines indicating the location of zero on each axis
#' @param nbin Binning resolution (you shouldn't have to mess with the default)
#' @param bandwidth Binning bandwidth (you shouldn't have to mess with the default)
#' @param cr Color ramp.  The only one supplied currently is the default
#' @param col Not sure about this one, but the default is ok
#' @param nrpoints Not sure about this one, but the default is ok
#' @param max_channel Max value the cytometer can produce.  See \code{\link{ax}}
#' @param tx tick style.  One of \code{'biexp', 'log', 'linear'}
#' @param ty tick style.  One of \code{'biexp', 'log', 'linear'}
#' @param plotaxt Logical, should we draw the axes?
#' @param ticksize Weight of the tickmarks in relative coordinates.
#' @param ... Additional graphical parameters passed to \code{flowCore \link[flowCore]{flowFrame-class}}
#' @examples
#'
#' # get some example data
#' filename = system.file("extdata", "example1.fcs", package = "wadeTools")
#' ff = get_sample(filename)
#'
#' # a plot with linear coordinate scaling
#' pplot(ff, parameters = c("FSC-A", "SSC-A"), tx = 'linear', ty = 'linear')
#'
#' # another plot, with the default biexponential scaling
#' pplot(ff, c("CD4PETR", "CD8Q705"))
#' @export
pplot <- function (ff, parameters, blueBackground=FALSE, showZero=FALSE, nbin=501, bandwidth=0.001,
                   cr = blob_color(blueBackground=FALSE), col = "black", nrpoints = 0,
                   max_channel = 262143,
                   tx=c("biexp", "asinh", "log", "linear"), ty=c("biexp", "asinh", "log", "linear"),
                   plotaxt = TRUE,  ticksize=1, ...) {

  requireNamespace("flowCore")
  requireNamespace("fields")

  if(!is.null(cr)){
    suppressWarnings (flowCore::plot (ff, parameters, colramp=cr, nbin=nbin, band=bandwidth, nrpoints=nrpoints, axes=FALSE, ...))
  }
  else{
    #       suppressWarnings (plot (ff, parameters, col=col, nbin=nbin, band=bandwidth, nrpoints=nrpoints, axes=FALSE, ...))

    suppressWarnings(plot (exprs(ff)[,parameters[1]],exprs(ff)[,parameters[2]], pch=20, col=col,cex =.2,  axes=F,xlab=parameters[1],ylab=parameters[2], ...))
    box()
  }

  if(plotaxt==TRUE){
    suppressWarnings (ax (1, max_channel = max_channel, type=tx, ticksize=ticksize))
    suppressWarnings (ax (2, max_channel = max_channel, type=ty, ticksize=ticksize))
  }
  if (showZero) {
    xline(0, lty='dotdash')
    yline(0, lty='dotdash')
  }
}

# for now this assumes biexp transform if tx != linear (or ty)
pplot.with.hist = function (ff, parameters, tx='biexp', ty='biexp', mthresh=c(.001, .001), bandwidth=c(.03, .03), find.thresholds=FALSE) {

  requireNamespace("fields")
  requireNamespace("KernSmooth")

  # set up the layout
  laymat = matrix (c(1, 3, 0, 2), nrow=2, ncol=2)
  widths = c(1, .3)
  heights = c(.3, 1)

  layout (laymat, widths=widths, heights=heights)


  pnames = colnames(ff)
  pind1 = grep (parameters[1], pnames)
  pind2 = grep (parameters[2], pnames)
  kde1 = bkde (exprs(ff)[,pind1], bandwidth=bandwidth[1], gridsize=2001)
  kde2 = bkde (exprs(ff)[,pind2], bandwidth=bandwidth[2], gridsize=2001)
  mn1 = parameters(ff)$minRange[pind1]
  mn2 = parameters(ff)$minRange[pind2]
  mx1 = parameters(ff)$maxRange[pind1]
  mx2 = parameters(ff)$maxRange[pind2]

  if (find.thresholds) {
    loc.min.x = find.local.minima (kde1, thresh=mthresh[1])$x
    loc.min.y = find.local.minima (kde2, thresh=mthresh[2])$x

  }


  par(mar=c(0,4,0,0)+.1)
  npts = length(kde1$x)
  kde1$y[1] = kde1$y[npts] = 0
  plot (kde1, type='l', xlim=c(mn1, mx1), xaxt='n', yaxt='n', xlab='', ylab='')
  polygon (kde1$x, kde1$y, col='gray')
  if (find.thresholds) {
    xline(loc.min.x, lty='dotdash')
  }
  # need to rotate the y histogram
  par (mar=c(4, 0, 0, 0)+.1)
  npts = length(kde2$x)
  kde2$y[1] = kde2$y[npts] = 0
  plot (kde2$y, kde2$x, type='l', ylim=c(mn2, mx2), xaxt='n', yaxt='n', xlab='', ylab='')
  polygon (kde2$y, kde2$x, col='gray')
  if (find.thresholds) {
    yline(loc.min.y, lty='dotdash')
  }

  par(mar=c(4, 4, 0, 0)+.1)
  pplot (ff, parameters, tx=tx, ty=ty, xlim=c(mn1,mx1), ylim=c(mn2,mx2))
  if (find.thresholds) {
    xline(loc.min.x, lty='dotdash', lwd=3)
    yline(loc.min.y, lty='dotdash', lwd=3)
  }

  # reset par to reasonable values
  par (mfrow=c(1,1), mar=c(5,4,4,1))

  if (find.thresholds) {
    invisible (list(x.thresh=loc.min.x, y.thresh=loc.min.y))
  } else {
    invisible()
  }
}

pplot.index = function (ff, parameter, idx.or.time=c("idx", "time"), tx='linear', ty='biexp', ...) {
  idx.or.time = match.arg(idx.or.time)
  if (idx.or.time == "idx") {
    index=1:nrow(ff)
    # incorporate index as a new parameter in a tmp flowFrame
    tmpmat = exprs(ff)
    tmpmat = cbind (tmpmat, index=index)
    tmp = flowFrame (tmpmat)
    pplot (tmp, plist=c("index", parameter), tx=tx, ty=ty, ...)
  } else {
    tparam = colnames(ff)[which(tolower(colnames(ff)) == "time")]
    pplot (ff, plist=c(tparam, parameter), tx=tx, ty=ty, ...)
  }
}

pplot.show.bad = function (ff, parameter, idx.or.time=c("idx", "time"), tx='linear', ty='biexp', print.lost=TRUE, ...) {
  idx.or.time = match.arg(idx.or.time)
  pplot.index (ff = ff, parameter = parameter, idx.or.time = idx.or.time, tx=tx, ty=ty, ...)
  tot.events = nrow(ff)

  # find events with clean == 0
  if ("clean" %in% colnames(ff)) {
    tmp = Subset (ff, rectangleGate("clean"=c(-1e-1, 1e-1)))
    lost.events = nrow(tmp)

    xpos = .4
    if (idx.or.time == "idx") {
      xvals = which (exprs(ff)[,"clean"] == 0.0)
      points (xvals, exprs(tmp)[,parameter], pch=20, cex=.2, col='black')
      xtext = xpos * nrow(ff)
    } else {
      xp = colnames(ff)[which(tolower(colnames(ff)) == "time")]
      points (exprs(tmp)[,xp], exprs(tmp)[,parameter], pch=20, cex=.2, col='black')
      xtext = xpos * max(exprs(ff)[,xp])
    }
  } else {
    lost.events = 0
  }
  ytext = 0.8 * max (exprs(ff)[,parameter])
  if (print.lost) {
    text (xtext, ytext, labels = sprintf ("Eliminated %d events  (%.2f%%)", lost.events, 100*lost.events/tot.events), pos=4)
  }
}



