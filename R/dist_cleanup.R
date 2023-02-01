# cleanup.R
#
# Functions for cleaning up flowFrames.
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without Still Pond Cytomics express written consent.             ##
################################################################################
################################################################################
#



#' @title Get Rid of Events on the Rail
#' @description Derail (that is, remove events on the positive axis) of selected
#' parameters
#' @param ff The input flowFrame
#' @param parameters The parameters of ff you wish to derail
#' @param fac A factor < 1.0 to use as a 'buffer zone'
#' @details derail consults the **maxRange** values in the phenodata slot of the
#' flowFrame ff.  If something has gone wrong, it's likely due to the transformation
#' or normalization of data in the flowFrame.
#' @return A de-railed flowFrame.
#' @export
derail = function(ff, parameters=NULL, fac=0.99) {
  if (is.null(parameters)) {
    stop ("You must supply at least one parameter to derail")
  }

  cnames = colnames(exprs(ff))

  # parameters must all be in ff
  tmp = parameters %in% cnames
  if (sum(tmp) != length (tmp)) {
    stop ("Parameters must all be in ff")
  }
  ind = which (cnames %in% parameters)
  if (is (ff, "flowSet")) {
    mx = parameters(ff[[1]])$maxRange
  } else if (is(ff, "flowFrame")) {
    mx = parameters(ff)$maxRange
  } else {
    stop ("First arg must be either a flowframe or a flowSet")
  }
  mx = fac * mx
  myexpr = vector('character')
  k = 1
  for (i in ind) {
    myexpr[k] = paste ("'", cnames[i], "'=c(-Inf, ", mx[i], ")", sep="")
    k = k + 1
  }
  tmp = paste (myexpr, collapse=", ")
  gexpr = paste ("list (", tmp, ")", sep="")

  out = Subset (ff, rectangleGate (filterId="derail", .gate=eval(parse(text=gexpr))))

  out
}

#' @title Get Rid of Events Whose Parameter Value is Below Zero
#' @description Eliminate values less than min.value from a flowFrame.
#' @param ff The input flowFrame
#' @param parameters The parameters of ff you wish to derail
#' @param min.value The value below which events will be removed from the flowFrame
#' @return A de-neg'd flowFrame.
#' @export
deneg = function (ff, parameters=NULL, min.value=0) {
  if (is.null(parameters)) {
    stop ("You must supply at least one parameter to derail")
  }

  cnames = colnames(ff)

  # parameters must all be in ff
  tmp = parameters %in% cnames
  if (sum(tmp) != length (tmp)) {
    stop ("Parameters must all be in ff")
  }
  ind = which (cnames %in% parameters)

  myexpr = vector('character')

  k = 1
  for (i in ind) {
    myexpr[k] = paste ("'", cnames[i], "'=c(", min.value, ", Inf)", sep="")
    k = k + 1
  }
  tmp = paste (myexpr, collapse=", ")
  gexpr = paste ("list (", tmp, ")", sep="")

  out = Subset (ff, rectangleGate (filterId="deneg", .gate=eval(parse(text=gexpr))))

  out
}

