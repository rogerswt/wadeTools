# blob_utils.R
#
# Blob utility functions.
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#

# draw a circle around the center of a blob
circle <- function (blob, radius, length=100) {
  cen1 <- mean(blob[,1])
  cen2 <- mean(blob[,2])

  out <- matrix(nrow=(length+1), ncol=2)
  colnames(out) <- colnames(blob)

  for (i in 1:length) {
    theta <- 2*pi*i/length
    out[i,1] <- cen1 + radius * cos(theta)
    out[i,2] <- cen2 + radius * sin(theta)
  }
  out[length+1,] <- out[1,]	# join the beginning and end (aesthetics)
  out
}

# draw a circle around a point
circle.point = function(pt, radius, length=100) {
  cen1 <- pt[1]
  cen2 <- pt[2]

  out <- matrix(nrow=(length+1), ncol=2)
  colnames(out) <- names(pt)

  for (i in 1:length) {
    theta <- 2*pi*i/length
    out[i,1] <- cen1 + radius * cos(theta)
    out[i,2] <- cen2 + radius * sin(theta)
  }
  out[length+1,] <- out[1,]	# join the beginning and end (aesthetics)
  out
}

# helper function
close.contour = function (blob) {
  x = blob[,1]
  y = blob[,2]
  npts = nrow(blob)
  if ((x[1] != x[npts]) | (y[1] != y[npts])) {
    blob = rbind(blob, c(x[1], y[1]))
  }
  blob
}
#' @title Get the Convex Hull of a Blob
#' @description find the convex hull of a polygon
#' @param blob A polygon, for example the result of blob.boundary()
#' @return Another polygon
#' @export
get.hull <- function (blob) {
  # 2018-08-27 WTR - handle edge case that blob is a single point
  if (is.vector(blob)) {
    cnames= names(blob)
    blob = matrix(blob, nrow = 1, ncol = 2)
    colnames(blob) = cnames
    return(blob)
  }
  if(!is.matrix(blob)) {
    cnames = names(blob)
    blob = matrix(blob, ncol = 2)
    colnames(blob) = cnames
    return(blob)
  }
  x <- blob[,1]
  y <- blob[,2]
  hull <- chull(x, y)
  poly <- matrix(c(x[hull], y[hull]), nrow=length(hull), ncol=2)
  # get rid of dupes
  poly <- unique(poly)
  # close the contour
  poly <- rbind(poly, poly[1,])
  colnames(poly) <- colnames(blob)
  return(poly)
}

#' @title Inflate a Blob Contour
#' @description Expand the input blob radially by a constant distance.
#' @param blob A polygon, for example the result of blob.boundary()
#' @param dist The distance by which to expand the blob
#' @return Another polygon
#' @details The input blob should be convex.  You can use get.hull() first if you're not sure.
#' @export
inflate.contour <- function (blob, dist) {
  require("splancs")
  x <- blob[,1]
  y <- blob[,2]
  len <- length(x)
  x_infl <- mat.or.vec (len, 1)
  y_infl <- mat.or.vec (len, 1)

  for (i in 1:len - 1) {
    dx <- x[i+1] - x[i]
    dy <- y[i+1] - y[i]
    ux <- -dy / sqrt(dx*dx + dy*dy)
    uy <-  dx / sqrt(dx*dx + dy*dy)
    x_infl[i] <- x[i] + ux*dist
    y_infl[i] <- y[i] + uy*dist
  }
  # do the closure point
  x_infl[len] <- x_infl[1]
  y_infl[len] <- y_infl[1]

  new_poly <- as.points (x_infl, y_infl)
  colnames(new_poly) <- colnames(blob)
  return (new_poly)

}

#' @title Smooth a Blob Contour
#' @description Create a smoother representation of the input blob.
#' @param blob A polygon, for example the result of blob.boundary()
#' @param npts The number of points used in smoothing.  The larger the number
#' the smoother the result.
#' @return Another polygon
#' @export
smooth.contour <- function (blob, npts=5) {
  require("splancs")
  x <- blob[,1]
  y <- blob[,2]
  len <- length(x)
  x_smo <- mat.or.vec (len, 1)
  y_smo <- mat.or.vec (len, 1)

  # prepare for circular wrapping
  xe <- mat.or.vec (len + npts - 1, 1)
  ye <- mat.or.vec (len + npts - 1, 1)
  b <- floor(npts/2)

  xe[(b+1):(len+b)] <- x
  ye[(b+1):(len+b)] <- y
  for (i in 1:(b)) {
    xe[i] <- x[len - i]
    ye[i] <- y[len - i]
    xe[len+b+i] <- x[i]
    ye[len+b+i] <- y[i]
  }

  # now do the smoothing
  for (i in 1:len) {
    j <- i+b
    x_smo[i] <- sum(xe[(j-b):(j+b)]) / npts
    y_smo[i] <- sum(ye[(j-b):(j+b)]) / npts
  }

  new_poly <- as.points (x_smo, y_smo)
  colnames(new_poly) <- colnames(blob)
  return (new_poly)

}

# Worker function - what is the centroid of a contour? Assumes cont is a matrix
centroid = function (cont) {
  nr = nrow(cont)
  if ((cont[1,1] == cont[nr,1]) && (cont[1,2] == cont[nr,2])) {
    cont = cont[1:(nr-1),]
  }

  vals = c(mean(cont[,1], na.rm=T), mean(cont[,2], na.rm=T))
  names(vals) = colnames(cont)
  return (vals)
}

# Worker function - is a point inside a contour?
inside = function (p, contour) {
  sum_theta = 0
  for (i in 1:(nrow(contour)-1)) {
    x1 = contour[i,1] - p[1]
    y1 = contour[i,2] - p[2]
    x2 = contour[i+1,1] - p[1]
    y2 = contour[i+1,2] - p[2]
    sum_theta = sum_theta + angle2D ( x1, y1, x2, y2)
  }

  sum_theta = abs(sum_theta)
  crit = .1
  if (abs(sum_theta - 2*pi) < crit) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}

# worker function - what is the angle between two vectors?
angle2D = function (x1, y1, x2, y2) {
  theta1 = atan2 (y1, x1)
  theta2 = atan2 (y2, x2)
  dtheta = theta2 - theta1
  if (dtheta > pi) {dtheta = dtheta - 2*pi}
  if (dtheta < -pi) {dtheta = dtheta + 2*pi}

  dtheta
}

# helper function
cont2mat = function (cont, param) {
  nrows <- length(cont$x)
  mat <- matrix(NA, nrow=nrows, ncol=2)
  mat[,1] <- cont$x
  mat[,2] <- cont$y
  if ((mat[1,1] != mat[nrows,1]) && (mat[1,2] != mat[nrows,2])) {
    mat = rbind (mat, c(cont$x[1], cont$y[1]))  # close the contour
  }
  colnames(mat) <- param
  mat
}

#' @title Calculate a Blob's Area
#' @description Calculate the area of a blob (polygon).
#' @param blob A polygon, for example the result of blob.boundary()
#' @return The area of the polygon.
#' @details See http://stackoverflow.com/questions/16285134/calculating-polygon-area.
#' @export
polygon.area = function(blob) {
  if (!is.matrix(blob)) {blob = cont2mat(blob, param=c("x", "y"))}

  area = 0.0
  npts = nrow(blob)

  j = npts
  for (i in 1:npts) {
    area = area + (blob[j,1] + blob[i,1]) * (blob[j,2] - blob[i,2])
    j = i
  }
  area = area / 2.0

  area
}

# convexity defined as the ratio of polygon area to the area of its convex hull
#' @title Calculate a Blob's Convexity
#' @description Calculate the convexity of a blob (polygon).
#' @param blob A polygon, for example the result of blob.boundary()
#' @return The convexity of the polygon.
#' @details Convexity is defined as the ratio of the area of a polygon to the area
#' of its convex hull.  A convex polygon thus has a convexity of 1.0.
#' @export
polygon.convexity = function(blob) {
  a1 = polygon.area(blob)
  a2 = polygon.area(get.hull(blob))

  convexity = a1 / a2

  convexity
}
