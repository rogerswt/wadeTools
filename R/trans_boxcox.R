#
#  Box-Cox Transform
#  Lo et. al., Cytometry A 73A: 321-332, 2008
#


boxcox.transform  <- function (x, lambda, scale=NA, jitter=TRUE) {

  # add a small amount of random noise to smooth things out at low vals
  if (jitter) {
    noise <- rnorm (length(x), sd=1)
    x <- x + noise
  }

  if (is.na(scale)) {
    scale = 5.4 / ((262143^lambda - 1)/lambda)
  }

  val = scale * sign(x) * (abs(x)^lambda - 1) / lambda

  val
}


############################################################
#  This function should be exposed - flowCore compatibility
############################################################
#' @title Create a FlowCore Compatible Transform Object
#' @description Do the boxcox mathematical transform
#' @param transformId FlowCore compatible transformId
#' @param lambda boxcox lambda parameter
#' @param scale boxcox scale parameter
#' @param jitter Logical, should small values be jittered for visual appeal?
#' @description This function creates a flowCore compatible transform object.
#' @return an object of type \code{transform}, to be used with flowCore functions.
#' @export
boxcoxTransform <- function(transformId = "myboxcox", lambda = 0.1, scale = NA, jitter = FALSE) {
  t <- new("transform", .Data = function (x) boxcox.transform(x, lambda, scale, jitter))
  t@transformationId = transformId
  t
}

