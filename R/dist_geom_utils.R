#
# geom_utils.R
#
# Rotation of 2-d coordinates to help blob.boundary
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#


# rotate 2D coordinates about a center point
# mat is a 2-column matrix of coordinates
rotate = function(mat, center = c(0, 0), degrees = 0) {
  if (!is.matrix(mat)) {
    stop("rotate: mat must be a matrix")
  }
  if (ncol(mat) != 2) {
    stop("rotate: mat must have 2 columns")
  }

  x = mat[, 1]
  y = mat[, 2]
  xc = center[1]
  yc = center[2]
  r = sqrt((x - xc)^2 + (y - yc)^2)
  theta = atan2(y - yc, x - xc)
  xprime = (r * cos( theta + degrees * pi / 180)) + xc
  yprime = (r * sin(theta + degrees * pi / 180)) + yc

  out = matrix(c(xprime, yprime), ncol = 2)
  colnames(out) = colnames(mat)

  out
}

