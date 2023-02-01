#
# gaussian_utils.R
#
# gaussian utility functions.
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without Still Pond Cytomics express written consent.             ##
################################################################################
################################################################################
#

#' @title Find the Crossover Point of Two Gaussians
#' @description Find the crossover of 2 gaussians, that is, where they intersect.
#' @param l1 Lambda of the first gaussian
#' @param s1 Sigma of the first gaussian
#' @param m1 Mu of the first gaussian
#' @param l2 Lambda of the second gaussian
#' @param s2 Sigma of the second gaussian
#' @param m2 Mu of the second gaussian
#' @return The x-coordinate of the intersection.
#' @export
gaussian_crossover = function(l1, s1, m1, l2, s2, m2) {
  k = 2 * log((l1 * s2) / (l2 * s1))
  if (s1 == s2) {
    x = (-k * s2 ^ 2 + m1 ^ 2 - m2 ^ 2) / (2 * (m1 - m2))
  } else {
    x1 = (sqrt(-1 * s1 ^ 2 * s2 ^ 2 * (k * s1 ^ 2 - k * s2 ^ 2 - m1 ^ 2 + 2 * m1 * m2 - m2 ^ 2)) - m1 * s2 ^ 2 + m2 * s1 ^ 2) / (s1 ^ 2 - s2 ^ 2)
    x2 = (-sqrt(-1 * s1 ^ 2 * s2 ^ 2 * (k * s1 ^ 2 - k * s2 ^ 2 - m1 ^ 2 + 2 * m1 * m2 - m2 ^ 2)) - m1 * s2 ^ 2 + m2 * s1 ^ 2) / (s1 ^ 2 - s2 ^ 2)
    # which root gives intersection between the two means?
    t1 = ifelse(x1 > m1 & x1 < m2, TRUE, FALSE)
    t2 = ifelse(x2 > m1 & x2 < m2, TRUE, FALSE)
    if (t1 & !t2) {x = x1}
    if (t2 & !t1) {x = x2}
    if (t1 & t2)  {
      x = max(x1, x2)
      cat("both roots between gaussians\n")
    }
    if (!t1 & !t2) {
      x = bx(800)
      cat("neither root between gaussians\n")
    }
  }
  x
}

#' @title Generate a Gaussian From a KDE
#' @description Generate a gaussian curve, in the same form as our familiar
#' kde's (a 'gde'), from an input kde.  The resulting gaussian  intercepts the
#' input kde at the given height.
#' @param kde The input kde
#' @param height The height at which the generated gde
#' @param show Logical.  Should we show our work?
#' @return A gde.
#' @export
fit_gaussian_base = function(kde, height, show=FALSE) {
  kde = normalize.kde(kde)
  res = find.width(kde, frac = height)
  base = res$right - res$left
  gwidth = interpolate_width_gaussian(height)
  sigma = base / gwidth
  mu = res$center
  gde = my.gaussian(x = kde$x, mu = mu, sigma = sigma, scale = 'unit.height')

  if (show) {
    plot(kde, xlim = c(-.5, bx(1000)), xaxt = 'n')
    bax()
    lines(gde, col = 'red')
  }
  gde
}

#' @title Get the Gaussian Parameters of a GDE
#' @description Find sigma, mu for a given Gaussian Density Estimate
#' (see for example \code{\link{fit_gaussian_base}})
#' @param gde The input Gaussian Density Estimate
#' @return A list with 2 elements:
#' describe{
#'    \item{sigma}{The sigma value for the input gde}
#'    \item{mu}{ The mu value for the input gde}
#' }
#' @export
retrieve_sigma_mu = function(gde) {
  height = 0.5
  gde = normalize.kde(gde)
  res = find.width(gde, frac = height)
  base = res$right - res$left
  gwidth = interpolate_width_gaussian(height)
  sigma = base / gwidth
  mu = res$center
  return(list(sigma = sigma, mu = mu))
}

cumm_prob_thresh = function(gde, cumm.prob=0.9999) {
  cumm = cumm.kde(gde)
  idx = max(which(cumm$y < cumm.prob))
  thresh = cumm$x[idx]

  thresh
}

# tabulate the width of a gaussian at various heights.
# This width is relative to sigma.  In other words, the generated gaussian
# has sigma = 1.
tabulate_width_gaussian = function(height = seq(.001, .999, length = 999)) {
  kde = my.gaussian(x = seq(-10, 10, length = 6001), mu = 0, sigma = 1, scale = 'unit.height')
  width = vector(mode = 'numeric')
  k = 1
  for (h in height) {
    res = find.width(kde, frac = h)
    width[k] = res$right - res$left
    k = k + 1
  }
  return(list(h = height, w = width))
}

# tab is the result of tabulate_width_gaussian
interpolate_width_gaussian = function(height) {
  tab = tabulate_width_gaussian()
  idx = which(abs(tab$h - height) == min(abs(tab$h - height)))

  tab$w[idx]
}

# generate a gaussian
#' @export
my.gaussian = function(x, mu, sigma, scale=c("unit.height", "unit.area")) {
  g = exp(-.5 * ((x - mu) / sigma) ^ 2)

  if (scale == 'unit.area') {
    g = g / (sigma * sqrt(2 * pi))
  }

  return(list(x = x, y = g))
}

# Generate a gaussian kde function
#    mu = mean
#    sigma = sigma
#    lambda = integral
#' @title Generate a Gaussian GDE Distribution
#' @description Generate a Gaussian Density Function (gde) given mu, sigma and lambda.
#' @param lambda Lambda
#' @param sigma Sigma
#' @param mu Mu
#' @param x The x values at which the gaussian densities will be generated
#' @return A list with 2 elements:
#' \describe{
#'    \item{x}{The x values of the gde}
#'    \item{y}{The y values of the gde}
#' }
#' @export
gaussian_dist = function(lambda, sigma, mu, x = seq(-1, 5, length.out = 200)) {
  a = lambda / (sigma * sqrt(2 * pi))
  e = -0.5 * ((x - mu) / sigma) ^ 2
  g = a * exp(e)

  return(list(x = x, y = g))
}

bax = function(axis = 1) {
  ax(axis = axis, type = 'biexp')
}

