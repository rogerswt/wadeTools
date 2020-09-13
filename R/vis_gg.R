# vis_gg.R
#
#	The function ggflow() uses ggplot to display bivariate flow cytometry data.  It's a more
#	sophisticated and flexible replacement for pplot().
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomics LLC 2020.                 ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#
#
#' @title Draw a FlowJo-style (sort of) Plot Using ggplot Semantics
#' @description Draw a picture, with default colors reminiscent of FlowJo.
#' @param ff The flowFrame to be plotted
#' @param params A vector of length 2 indicating the parameters of ff to be plotted
#' @param colors A color palette for rendering densities of events.  Currently, one
#' of c("flowjo", "viridis", "plasma", "magma").  Can be abbreviated.
#' @param resolution One of c("medium", "coarse", "fine").  Can be abbreviated.
#' @param trans_sc The transformation that was applied to the scattering parameters of ff.
#' Currently, one of c("linear", "biexp", "asinh", "log").
#' @param trans_fl The transformation that was applied to the fluorescence parameters of ff.
#' Currently, one of c("linear", "biexp", "asinh", "log").
#' @param xlim The limits of the plot in the x direction. NULL will apply sensible defaults.
#' @param ylim The limits of the plot in the y direction. NULL will apply sensible defaults.
#' @return A ggplot object.
#' @description ggflow uses ggplot semantics to render bivariate flow cytometry plots.
#' If called "bare", it draws a picture.  If its return value is captured in a variable,
#' further manipulation of the plot is possible.  See the examples below.
#'
#' @examples
#'
#' # get some example data
#' filename = system.file("extdata", "example1.fcs", package = "wadeTools")
#' ff = get_sample(filename)
#' # Recall that get_sample by default applies linear transformation to scattering
#' # parameters and biexponential transformation to fluorescence parameters.
#'
#' # Plot ff on default FSC-A, SSC-A dimensions (notice linear axes)
#' ggflow(ff)
#'
#' # another plot of fluorescence parameters.  Note the biexponential axes.
#' ggflow(ff, c("CD4PETR", "CD8Q705"))
#'
#' # Add a title
#' ggflow(ff, c("CD4PETR", "CD8Q705")) + labs(title = "Title of Plot")
#'
#' # or if you prefer,
#' p = ggflow(ff, c("CD4PETR", "CD8Q705"))
#' p + labs(title = "Another plot with a title")
#'
#' # Finally, some of the other color schemes
#' library("gridExtra")
#' p1 = ggflow(ff, c("SSC-A", "LIVEDEAD")) + labs(title = "flowjo") + theme(legend.position = 'right')
#' p2 = ggflow(ff, c("SSC-A", "LIVEDEAD"), colors = "v") + labs(title = "viridis") + theme(legend.position = 'right')
#' p3 = ggflow(ff, c("SSC-A", "LIVEDEAD"), colors = "m") + labs(title = "magma") + theme(legend.position = 'right')
#' p4 = ggflow(ff, c("SSC-A", "LIVEDEAD"), colors = "p") + labs(title = "plasma") + theme(legend.position = 'right')
#' p = grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
#' p
#'
#' @export
ggflow = function(ff,
                  params = c("FSC-A", "SSC-A"),
                  colors = c("flowjo", "viridis", "plasma", "magma"),
                  resolution = c("medium", "coarse", "fine"),
                  trans_sc = c("linear", "biexp", "asinh", "log"),
                  trans_fl = c("biexp", "asinh", "log", "linear"),
                  xlim = NULL, ylim = NULL) {
  requireNamespace("ggplot2")
  requireNamespace("ggcyto")
  requireNamespace("viridis")
  colors = match.arg(colors)
  trans_sc = match.arg(trans_sc)
  trans_fl = match.arg(trans_fl)
  resolution = match.arg(resolution)

  if (resolution == "coarse") {
    bins = 500
  } else if (resolution == "medium") {
    bins = 750
  } else if (resolution == "fine") {
    bins = 1000
  }

  params_sc = detect_scatter_parameters(ff)
  params_fl = detect_fl_parameters(ff)

  # this crazy ".data" indexing is described in
  # https://dplyr.tidyverse.org/articles/programming.html
  p = ggplot(ff, mapping = aes(x = .data[[params[1]]], y = .data[[params[2]]]))
  binning = geom_bin2d(bins = bins, na.rm = TRUE)
  p = p + binning

  # apply a theme
  My_Theme = theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 24, hjust = 0.5),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 0),
    legend.position = 'none')
  p = p + My_Theme

  # construct appropriate axes
  params_sc = detect_scatter_parameters(ff)
  params_fl = detect_fl_parameters(ff)

  # x axis
  if (params[1] %in% colnames(ff)[params_sc]) {
    method = trans_sc
  } else {
    method = trans_fl
  }
  a = ticks_breaks_labels(ff, params[1], method = method)
  if (is.null(xlim)) {
    limits = a$range
  } else {
    limits = xlim
  }
  p = p + scale_x_continuous(breaks = a$major, limits = limits, minor_breaks = a$ticks, labels = a$labels)


  # y axis
  if (params[2] %in% colnames(ff)[params_sc]) {
    method = trans_sc
  } else {
    method = trans_fl
  }
  a = ticks_breaks_labels(ff, params[2], method = method)
  if (is.null(ylim)) {
    limits = a$range
  } else {
    limits = ylim
  }
  p = p + scale_y_continuous(breaks = a$major, limits = limits, minor_breaks = a$ticks, labels = a$labels)


  # Choose a color mapping
  if (colors == "flowjo") {
    p = p + scale_fill_gradientn(colours = hsv(h = seq(.6667, 0, length.out = 11)))
  }
  else {
    p = p + viridis::scale_fill_viridis(option = colors)
  }

  return(p)
}

ticks_breaks_labels = function(ff, param, method = c("biexp", "asinh", "log", "linear")) {
  method = match.arg(method)
  start_decade = 2
  decades = 8
  major <- (10^(start_decade:decades))
  all.ticks <- NULL
  for (i in 1:(length(major) - 1)) {
    all.ticks <- c(all.ticks, seq(major[i], major[i + 1], l = 10))
  }
  if (method == "biexp") {
    major = bx(major)
    all.ticks <- bx(all.ticks)
    labels = tick.labels(start_decade, decades)
    range = c(0, bx(2^18))
  } else if (method == "asinh") {
    major = w.arcsinh(major)
    all.ticks <- w.arcsinh(all.ticks)
    labels = tick.labels(start_decade, decades)
    range = c(0, w.arcsinh(2^18))
  } else if (method == "log") {
    major = log10(major)
    all.ticks <- log10(all.ticks)
    labels = tick.labels(start_decade, decades)
    range = c(0, log10(2^18))
  } else if (method == 'linear') {
    idx = which(exprs(ff)[, param] >= 0)
    major = pretty(exprs(ff)[idx, param])
    ticks = major
    labels = major
    idx = which(param == colnames(ff))
    range = c(0, parameters(ff)$maxRange[idx])
  }
  return(list(major = major, ticks = all.ticks, labels = labels, range = range))
}
