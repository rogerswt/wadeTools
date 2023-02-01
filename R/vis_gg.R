# vis_gg.R
#
#	The function ggflow() uses ggplot to display bivariate flow cytometry data.  It's a more
#	sophisticated and flexible replacement for pplot().
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomics LLC 2020.                 ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without Still Pond Cytomics express written consent.             ##
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
#' @param indicate_zero Boolean, should we indicate the location of 0?
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
#' # one fluorescence (biexp) versus one scattering (linear)
#' ggflow(ff, c("SSC-A", "CD3Q605"))
#'
#' # Add a title
#' ggflow(ff, c("CD4PETR", "CD8Q705")) + labs(title = "Title of Plot")
#'
#' # or if you prefer,
#' p = ggflow(ff, c("CD4PETR", "CD8Q705"))
#' p + labs(title = "Another plot with a title")
#'
#' # Check out the other color schemes
#' library("gridExtra")
#' p1 = ggflow(ff, c("CD3Q605", "LIVEDEAD")) + labs(title = "flowjo") + theme(legend.position = 'right')
#' p2 = ggflow(ff, c("CD3Q605", "LIVEDEAD"), colors = "v") + labs(title = "viridis") + theme(legend.position = 'right')
#' p3 = ggflow(ff, c("CD3Q605", "LIVEDEAD"), colors = "m") + labs(title = "magma") + theme(legend.position = 'right')
#' p4 = ggflow(ff, c("CD3Q605", "LIVEDEAD"), colors = "p") + labs(title = "plasma") + theme(legend.position = 'right')
#' p = grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
#' p
#'
#' # Override plot limits and default resolution
#' ggflow(ff, c("CD45RAQ655","CD11BAPCCY7"), xlim = bx(c(-1e3,1e6)), ylim = bx(c(-1e3,1e6)), res = 'f')
#'
#' # make a plot and add some blob boundaries
#' params = c("SSC-A", "CD3Q605")
#' x_range = 0.5
#' p1 = ggflow(ff, params)
#' bb1 = blob.boundary(ff, c("SSC-A", "CD3Q605"), x_range = c(x_range, Inf), location = c(1.5, bx(8000)), height = .25)
#' b1 = geom_path(bb1, mapping = aes(x = `SSC-A`, y = CD3Q605, group = 1))
#' bb2 = blob.boundary(ff, c("SSC-A", "CD3Q605"), x_range = c(x_range, Inf), location = c(3, bx(500)), height = .25)
#' b2 = geom_path(bb2, mapping = aes(x = `SSC-A`, y = CD3Q605, group = 1))
#' bb3 = blob.boundary(ff, c("SSC-A", "CD3Q605"), x_range = c(x_range, Inf), location = c(0, 0), height = .25)
#' b3 = geom_path(bb3, mapping = aes(x = `SSC-A`, y = CD3Q605, group = 1))
#' vline = geom_vline(aes(xintercept = x_range), linetype = "dotted")
#' p1 + b1 + b2 + b3 + vline
#'
#' @export
ggflow = function(ff,
                  params = c("FSC-A", "SSC-A"),
                  colors = c("flowjo", "viridis", "plasma", "magma"),
                  resolution = c("medium", "very_coarse","coarse", "fine"),
                  trans_sc = c("linear", "biexp", "asinh", "log"),
                  trans_fl = c("biexp", "asinh", "log", "linear"),
                  indicate_zero = TRUE,
                  xlim = NULL, ylim = NULL) {
  requireNamespace("ggplot2")
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
  } else if (resolution == "very_coarse") {
    bins = 250
  }

  params_sc = detect_scatter_parameters(ff)
  params_fl = detect_fl_parameters(ff)

  # this crazy ".data" indexing is described in
  # https://dplyr.tidyverse.org/articles/programming.html
  dat = data.frame(exprs(ff), check.names = FALSE)
  p = ggplot(dat, mapping = aes(x = .data[[params[1]]], y = .data[[params[2]]]))
  binning = geom_bin2d(bins = bins, na.rm = TRUE)
  p = p + binning

  if (indicate_zero) {
    p = p + geom_hline(yintercept = 0, linetype = 'dotdash')
    p = p + geom_vline(xintercept = 0, linetype = 'dotdash')
  }

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

#' @title Support for ggplot-style axis labeling
#' @description Placing tick marks and labels appropriate for several popular transformations.
#' @param ff The flowFrame to be plotted
#' @param param The parameter of ff to be plotted
#' @param method The transformation that was applied to the scattering parameters of ff.
#' Currently, one of c("linear", "biexp", "asinh", "log", "linear").
#' @return A list comprised of:
#'            major = vector of locations of major ticks
#'            ticks = vector of locations of minor ticks
#'            labels = vector of major tick labels
#'            range = vector of length 2, indicating min and max of range.  Default values
#'            can be overridden if desired.
#' @description This function provides support for ggplot style axis generation.
#' It is used in ggflow(), but sometimes you want to make custom plots that respect
#' the transformation you've used on your data.
#'
#' @examples
#'
#' # get some example data
#' filename = system.file("extdata", "example1.fcs", package = "wadeTools")
#' ff = get_sample(filename)
#' # Recall that get_sample by default applies linear transformation to scattering
#' # parameters and biexponential transformation to fluorescence parameters.
#'
#' a = ticks_breaks_labels(ff, param = "CD3Q605")
#'
#' # make a plot of a kernel density estimate of a univariate parameter
#' kde_3 = normalize.kde(bkde(exprs(ff)[,"CD3Q605"], band = 0.1, grid = 1001))
#'
#' # limit range of kde for plotting
#' tmp = data.frame(kde_3)
#' tmp = tmp[tmp$x > a$range[1] & tmp$x < a$range[2], ]
#' p = ggplot(tmp, aes(x = x, y = y)) +
#'          geom_path() + xlab("") + ylab("") +
#'          labs(title = "CD3") +
#'          theme(plot.title = element_text(size = 30, hjust = 0.5))
#' xax = scale_x_continuous(breaks = a$major, limits = a$range, minor_breaks = a$ticks, labels = a$labels)
#' p + xax
#'
#'
#' @export
ticks_breaks_labels = function(ff, param, method = c("biexp", "asinh", "log", "linear")) {
  method = match.arg(method)

  neg_major = -(10^(2:10))
  pos_major = 10^(2:10)
  major = c(neg_major, pos_major)
  # all.ticks <- NULL
  # for (i in 1:(length(major) - 1)) {
  #   all.ticks <- c(all.ticks, seq(major[i], major[i + 1], l = 10))
  # }
  neg.ticks = vector(mode = 'numeric')
  for (i in 1:(length(neg_major) - 1)) {
    neg.ticks = c(neg.ticks, seq(neg_major[i], neg_major[i + 1], l = 10))
  }
  pos.ticks = vector(mode = 'numeric')
  for (i in 1:(length(pos_major) - 1)) {
    pos.ticks = c(pos.ticks, seq(pos_major[i], pos_major[i + 1], l = 10))
  }
  all.ticks = sort(unique(c(neg.ticks, pos.ticks)))
  if (method == "biexp") {
    labels = gg_tick_labels(major)
    major = bx(major)
    all.ticks <- bx(all.ticks)
    range = c(bx(-200), bx(2^18))
  } else if (method == "asinh") {
    neg_major = -(10^(0:10))
    pos_major = 10^(0:10)
    major = c(neg_major, pos_major)
    neg.ticks = vector(mode = 'numeric')
    for (i in 1:(length(neg_major) - 1)) {
      neg.ticks = c(neg.ticks, seq(neg_major[i], neg_major[i + 1], l = 10))
    }
    pos.ticks = vector(mode = 'numeric')
    for (i in 1:(length(pos_major) - 1)) {
      pos.ticks = c(pos.ticks, seq(pos_major[i], pos_major[i + 1], l = 10))
    }
    all.ticks = sort(unique(c(neg.ticks, pos.ticks)))
    labels = gg_tick_labels(major)
    major = asx(major)
    all.ticks <- asx(all.ticks)
    range = c(asx(0), asx(2^14))
  } else if (method == "log") {
    labels = gg_tick_labels(major)
    major = log10(major)
    all.ticks <- log10(all.ticks)
    range = c(log10(1), log10(2^18))
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
