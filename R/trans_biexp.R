#
# trans_biexp.R
#
# Transform values in a flowFrame object for "biexponential" plotting
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#

############################################################
#	The following functions are utilities
#	and should NOT be exposed
############################################################

kernel_function <- function(x, a = 0.002) {
    return (log(0.5*(x + sqrt((x)^2 +1/a^2)))/log(10))
}

biexp.transform  <- function(x, a = 0.002, full_scale = 262143, jitter = TRUE) {

  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
    # add a small amount of random noise to smooth things out at low vals
	if (jitter) {
		noise <- rnorm(length(x), sd = 1)
		x <- x + noise
	}

	tmp <- kernel_function(x, a)

	# offset downwards, then compensate by multiplying asymptotic value
	offset <- kernel_function(0, a)
	intercept = log10(full_scale)
	fac <- intercept / (kernel_function(full_scale, a) - offset)

	tmp <-  fac * (tmp - offset)
	return(tmp)
}


inv_kernel_function <- function(x, a = 0.002) {
	out <- .25 * 10 ^ (-x) * (-1 / a ^ 2 + 4 * 10 ^ (2 * x))
	out
}

inv.biexp.transform <- function(x, a = .002, full_scale = 262143) {

	# reverse the maneuvers
	offset <- kernel_function(0, a)
	intercept = log10(full_scale)
	fac <- intercept / (kernel_function(full_scale, a) - offset)

	tmp <- inv_kernel_function(x/fac + offset, a)
	tmp
}

############################################################
#	This function should be exposed
############################################################
#' @title Calculate a Biexponential Transformation
#' @description Do the biexponential mathematical transform
#' @param ff a flowFrame
#' @param a The "a" parameter.  Do not override unless you know what you're doing.
#' @param params The parameters of ff to transform
#' @description This function performs the biexponential transform directly on ff.
#' See \code{\link{biexpTransform}} for flowCore compatibility.
#' @return A transformed flowFrame
#' @export
biexp <- function (ff, a = 0.002, params) {

	# if using symbolic names for params, convert to numeric
	if (is.character(params)) {
		params <- which(colnames(ff) %in% params)
	}

   # make a copy of the input object
   fout <- ff

   # compute transformed object
   for (parm in params) {
      exprs(fout)[,parm] <- biexp.transform (exprs(ff)[,parm], a)
   }

   # set the minRange and maxRange values appropriately
   parameters(fout)$minRange[params] <- biexp.transform(parameters(fout)$minRange[params])
   parameters(fout)$maxRange[params] <- biexp.transform(parameters(fout)$maxRange[params])

   return(fout)
}


############################################################
#	This function should be exposed
############################################################
#' @title Calculate an Inverse Biexponential Transformation
#' @description Do the inverse biexponential mathematical transform
#' @param ff a flowFrame
#' @param a The "a" parameter.  Do not override unless you know what you're doing.
#' @param params The parameters of ff to transform
#' @description This function performs the inverse biexponential transform directly on ff.
#' @return A transformed flowFrame
#' @export
inv.biexp <- function(ff, a = 0.002, params) {
	# if using symbolic names for params, convert to numeric
	if (is.character(params)) {
		params <- which(colnames(ff) %in% params)
	}

   # make a copy of the input object
   fout <- ff

   # compute transformed object
   for (parm in params) {
      exprs(fout)[,parm] <- inv.biexp.transform(exprs(ff)[,parm], a)
   }

   # set the minRange and maxRange values appropriately
   parameters(fout)$minRange[params] <- inv.biexp.transform(parameters(fout)$minRange[params])
   parameters(fout)$maxRange[params] <- inv.biexp.transform(parameters(fout)$maxRange[params])

   return(fout)
}


############################################################
#	This function should be exposed - flowCore compatibility
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
biexpTransform <- function (transformId = "mybiexp", a = 0.002, full_scale = 262143, jitter = TRUE) {
	t <- new("transform", .Data = function (x) biexp.transform(x, a, full_scale, jitter))
	t@transformationId = transformId
	t
}

#' @title Shorthand Function for Biexponential Transformation
#' @description A shorthand function to perform the biexponential transform on numeric values
#' @param x A numeric scalar or vector to be transformed
#' @return Transformed values.
#' @seealso \code{\link{ibx}}
#' @export
bx = function(x) {
  idx = which(x == -Inf)
  res = biexp.transform(x, jitter = FALSE)
  if  (length(idx) > 0) {
    res[idx] = -Inf
  }
  res
}

#' @title Shorthand Function for Inverse Biexponential Transformation
#' @description A shorthand function to perform the inverse biexponential transform on numeric values
#' @param x A numeric scalar or vector to be transformed
#' @return Transformed values.
#' @seealso \code{\link{bx}}
#' @export
ibx = function(x) {
  res = inv.biexp.transform(x)

  res
}

