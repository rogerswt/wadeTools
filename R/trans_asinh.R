#
# trans_asinh.R
#
# Transform values in a flowFrame object using hyperbolic arcsin
#
#    NOTE: the form chosen here is compatible with default values in CellEngine.
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without Still Pond Cytomics express written consent.             ##
################################################################################
################################################################################
#

############################################################
#	The following functions are utilities
#	and should NOT be exposed
############################################################

asinh.transform = function(x, cofactor = 5, jitter = TRUE) {

  # add a small amount of random noise to smooth things out at low vals
  if (jitter) {
    noise <- rnorm(length(x), sd = 1)
    x <- x + noise
  }

  val = asinh(x/cofactor)

  val
}

inv.asinh.transform = function(x, cofactor = 5) {
  val = cofactor * sinh(x)

  val
}

############################################################
#	These functions should be exposed
############################################################

#' @title Calculate a Arcsinh Transformation
#' @description Do the arcsinh mathematical transform
#' @param ff a flowFrame
#' @param a The "a" parameter.  Do not override unless you know what you're doing.
#' @param params The parameters of ff to transform
#' @description This function performs the arcsinh transform directly on ff.
#' See \code{\link{biexpTransform}} for flowCore compatibility.
#' @return A transformed flowFrame
#' @export
w.arcsinh <- function(ff, cofactor = 5, params) {

  # if using symbolic names for params, convert to numeric
  if (is.character(params)) {
    params <- which(colnames(ff) %in% params)
  }

  # make a copy of the input object
  fout <- ff

  # compute transformed object
  for (parm in params) {
    exprs(fout)[,parm] <- asinh.transform(exprs(ff)[,parm], cofactor)
  }

  # set the minRange and maxRange values appropriately
  parameters(fout)$minRange[params] <- asinh.transform(parameters(fout)$minRange[params])
  parameters(fout)$maxRange[params] <- asinh.transform(parameters(fout)$maxRange[params])

  return(fout)
}

#' @title Calculate an Inverse Arcsinh Transformation
#' @description Do the inverse arcsinh mathematical transform
#' @param ff a flowFrame
#' @param a The "a" parameter.  Do not override unless you know what you're doing.
#' @param params The parameters of ff to transform
#' @description This function performs the inverse arcsinh transform directly on ff.
#' @return A transformed flowFrame
#' @export
w.inv.arcsinh <- function(ff, cofactor = 5, params) {
  # if using symbolic names for params, convert to numeric
  if (is.character(params)) {
    params <- which(colnames(ff) %in% params)
  }

  # make a copy of the input object
  fout <- ff

  # compute transformed object
  for (parm in params) {
    exprs(fout)[,parm] <- inv.arcsinh.transform(exprs(ff)[,parm], cofactor)
  }

  # set the minRange and maxRange values appropriately
  parameters(fout)$minRange[params] <- inv.arcsinh.transform(parameters(fout)$minRange[params])
  parameters(fout)$maxRange[params] <- inv.arcsinh.transform(parameters(fout)$maxRange[params])

  return(fout)
}

############################################################
#	flowCore compatibility
############################################################
#' @title Create a FlowCore Compatible Transform Object
#' @description Create a flowCore compatible transform object.
#' @param transformId FlowCore compatible transformId
#' @param a The "a" parameter.  Do not override unless you know what you're doing.
#' @param full_scale The value of full-scale events (default = 2^18 - 1)
#' @param jitter Logical, should small values be jittered for visual appeal?
#' @description This function creates a flowCore compatible transform object.
#' @return an object of type \code{transform}, to be used with flowCore functions.
#' @export
asinhTransform <- function (transformId = "myasinh", cofactor = 5, jitter = TRUE) {
  t <- new("transform", .Data = function (x) asinh.transform(x, cofactor, jitter))
  t@transformationId = transformId
  t
}

#' @title Shorthand Function for Arcsinh Transformation
#' @description A shorthand function to perform the arcsinh transform on numeric values
#' @param x A numeric scalar or vector to be transformed
#' @return Transformed values.
#' @seealso \code{\link{ibx}}
#' @export
asx = function(x, cofactor = 5) {
  idx = which(x == -Inf)
  res = asinh.transform(x, cofactor, jitter = FALSE)
  if  (length(idx) > 0) {
    res[idx] = -Inf
  }
  res
}

#' @title Shorthand Function for Inverse Arcsinh Transformation
#' @description A shorthand function to perform the inverse arcsinh transform on numeric values
#' @param x A numeric scalar or vector to be transformed
#' @return Transformed values.
#' @seealso \code{\link{bx}}
#' @export
iasx = function(x, cofactor = 5) {
  res = inv.asinh.transform(x, cofactor)

  res
}



