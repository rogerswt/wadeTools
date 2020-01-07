# blob_boundary.R
#
# gaussian utility functions.
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#


#' @title Find Blobs in 2 Dimensions
#' @description The idea is to compute the kernel density estimate of a 2D slice
#' of a flowFrame, compute contours at a specified height, and
#' then select that contour which is either (a) closest to a specified location,
#' or (b) is the 'strongest' peak in the flowFrame.  Optionally, you can specify
#' that the contour be convex, in which case the function will attempt to find
#' the largest convex contour at either location.
#' @param ff The input flowFrame
#' @param parameters Exactly 2 parameters included in ff
#' @param rotate Rotation angle (in degrees) of the blob you're looking for
#' @param location The target location of the blob (ignored if strongest = TRUE)
#' @param strongest Logical, whether to look for the blob containing the highest density
#' @param bandwidth 2D measure of bandwidth used to compute the Kernel Density Estimate
#' @param gridsize 2D measure of gridsize used to compute the Kernel Density Estimate
#' @param height The height at which to compute contours (height range is 0 - 1, 1 being the most dense peak)
#' @param convex Logical, whether to find the largest convex polygon at the target location
#' @param height1 Starting height for searching for convex blob (only used if convex = TRUE)
#' @param height2 Ending height for searching for convex blob (only used if convex = TRUE)
#' @param delta_h Step size used for searching for convex blob (only used if convex = TRUE)
#' @param log.transform Logical, should we use the density or its logarithm
#' @details blob.boundary uses \code{\link[grDevices]{contourLines}} to generate contours in density data.
#' It can get confused in several circumstances.  Sometimes it's
#' hard to detect a weak blob in the presence of a very strong blob, particularly
#' if they're close to each other and there's significant density connecting them.  It's
#' pretty easy to see why: contour lines will easily span between the two, effectively
#' connecting them into one bigger blob.  It is sometimes possible to pre-gate the data
#' to eliminate dense regions known to be uninteresting, and which may interfere
#' with blob finding.
#'
#' It is sometimes useful to compute contours on the **logarithm** of density, rather
#' than on the original density representation.  If you set log.transform = TRUE, it's
#' usually the case that a higher height parameter will work better than the smaller
#' (default) value.  This is because taking the log compresses the dynamic range, so
#' effectively everything is sort of 'high'.
#' @return A polygon, described as a 2-column matrix of coordinates.  You can use this polygon
#' to gate data using the **flowCore** function \code{polygonGate(.gate = bb)}.
#' @export
blob.boundary <- function (ff, parameters=c("FSC-A", "SSC-A"), rotate = 0.0,
                           location, strongest = FALSE,
                           bandwidth=c(.02, .02), gridsize=c(201, 201),
                           height=.1, convex = FALSE,
                           height1 = height, height2 = 0.05, delta_h = 0.01,
                           log.transform=FALSE) {
	require ("KernSmooth")
	require ("flowCore")

  if (!(is(ff)[[1]]) == "flowFrame") {
		stop ("first argument must be a flowFrame\n")
	}

	# extract a matrix of values
	mat <- exprs(ff)[,parameters]

	if (rotate != 0) {
	  mat = rotate(mat, center = location, degrees = rotate)
	}

	# compute a reasonable bandwith for the kde
	bw1 <- bandwidth[1] * max (mat[,1])
	bw2 <- bandwidth[2] * max (mat[,2])
	# do the kernel density estimate
	kde <- bkde2D (mat, bandwidth=c(bw1, bw2), gridsize=gridsize)

	# normalize the density estimate for sanity
	kde$fhat <- kde$fhat / max(kde$fhat)

  if (log.transform) {
    epsilon = 1e-4
    kde$fhat = log10(epsilon + kde$fhat)
    # renormalize to [0,1]
    mx = max(kde$fhat)
    mn = min(kde$fhat)
    kde$fhat = (kde$fhat - mn) / (mx - mn)
  }

	# compute contour lines at the given height
	cont <- contourLines (kde$x1, kde$x2, kde$fhat, levels=height)
	if (length(cont) == 0) {
		cat ("No contours were found!\n")
		return
	}

	if (strongest) {
	  # pick the contour that encloses the most events - 2018-08-12 WTR
	  nev = vector(mode = 'numeric')
	  for (i in 1:length(cont)) {
	    mat = cont2mat(cont[[i]], param = parameters)
	    nev[i] = nrow(Subset(ff, polygonGate(.gate = mat)))
	  }
	  idx = which(nev == max(nev))
	  out = cont2mat(cont[[idx]], param = parameters)
	} else {
	  # pick the  contour closest to the target location
	  out = nearest_contour(cont, location, parameters)
	}

	# the convex option searches for the largest contour that remains convex
	if (convex) {
	  tol = 1e-3   # tolerance for convexity
	  # update location based on what we've found
	  location = centroid(out)

	  # start high, search downwards
	  out = nearest_contour(contourLines(kde$x1, kde$x2, kde$fhat, levels = height1), location, parameters)
	  for (h in seq(height1, height2, by = -delta_h)) {
	    prev = out
	    cont <- contourLines(kde$x1, kde$x2, kde$fhat, levels = h)
	    out = nearest_contour(cont, location, parameters)
	    cx = polygon.convexity(out)
	    #DEBUG: cat("h =", h, ", convexity =", cx, "\n")
	    if (cx < 1.0 - tol) {
	      out = prev
	      break
	     }
	  }
	}

	if (rotate != 0) {
	  out = rotate(out, center = location, degrees = -rotate)
	}
	out
}

# helper function - find the contour whose center is nearest the specified location
# return in matrix form
nearest_contour = function(cont, location, parameters) {
  dist <- vector("numeric")
  for (i in 1:length(cont)) {
    xcen <- mean(cont[[i]]$x)
    ycen <- mean(cont[[i]]$y)
    dist[i] <- (xcen - location[1])^2 + (ycen - location[2])^2
  }
  nearest <- sort(dist, index.return=T)$ix[1]
  out = cont2mat(cont[[nearest]], param = parameters)

  out
}
