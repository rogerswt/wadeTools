#
# ff_importutils.R
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#


#' @title Read an FCS File and Do Some Processing
#' @description Read a FCS file using \code{\link[flowCore]{read.FCS}} and return a
#' flowFrame.  Along the way, optionally apply compensation, transformation,
#' removal of 'rail' events, and apply some parameter renaming for convenience
#' and nice pictures.
#' @param fn The fully qualified FCS filename
#' @param compensate Logical, apply the internal compensation matrix
#' @param transform Logical, apply linear transformation to scattering parameters
#' and biexponential transformation to fluorescence parameters
#' @param derail Logical, get rid of events on the FSC-A and SSC-A positive axes
#' @param nice.names Logical, swap 'desc' for 'name' in the flowFrame
#' @param verbose Logical, report progress on the console
#' @details This function is a convenience wrapper to apply commonly used processing
#' steps as FCS data are read in.  It's often important that all files in a project
#' get processed identically.  Consolidating most of the steps in one function helps
#' us avoid mistakes.
#' @return A flowFrame, properly processed.
#' @export
get_sample = function(fn, compensate=TRUE, transform=TRUE, derail=TRUE, nice.names = TRUE, verbose=FALSE) {


  requireNamespace("flowCore")

  ff = read.FCS(fn)
  if (verbose) {message("reading")}
  fl_params = which(flowCore::colnames(ff) %in% colnames(keyword(ff)$SPILL))
  sc_params = which(grepl(pattern = "FSC", x = flowCore::colnames(ff)) | grepl(pattern = "SSC", x = flowCore::colnames(ff)))

  if (compensate) {if (verbose) {message("compensating")}; ff = autocomp(ff)}

  if (derail) {
    if (verbose) {message("derailing")}
    ff = derail(ff = ff, parameters = c("FSC-A", "SSC-A"))
  }
  if (transform) {
    if (verbose) {message("transforming")}
    ff = doTransform(ff, cols = sc_params, method = 'linear')
    ff = doTransform(ff, cols = fl_params, method = 'biexp')
  }
  if (nice.names) {
    if (verbose) {message("nice names")}
    ff = swap.names(ff)
  }

  ff
}



# doTransform
#  This function performs biexponential transform on all parameters except SSC-W
# Jitter is set to false for biexponential to get reproducible results (otherwise would result in random #)
#  return value: a transformed flowFrame
#' @title A Convenience Wrapper for Transformation
#' @description Convenience wrapper for \code{\link[flowCore]{transformList}} and \code{\link[flowCore]{transform}}.
#' @param ff The flowFrame to be transformed
#' @param cols The columns (parameters) of the flowFrame to be transformed
#' @param method The transformation function. One of \code{biexp}, \code{log} or \code{linear}.
#' @param fac The scaling factor used only for \code{method = 'linear'}.
#' @details Please refer to the code for the math details.  The biexp method is specific to this
#' package.  It's similar to other biexponential transforms (e.g. flowCore, FlowJo) but differs
#' in detail.
#' @return The transformed flowFrame.
#' @export
doTransform <- function(ff,cols = c(1:5, 7:13),method = c("biexp","log","linear"), fac = 5.4/262143) {

  requireNamespace("flowCore")

  if (is(ff,"flowFrame")) {
    method = match.arg(method)
    if (method == "biexp") {
      bx <- biexpTransform(jitter = F)
      bxlist <- transformList (flowCore::colnames(ff)[cols], bx)
      return(flowCore::transform(ff, bxlist))
    }
    if (method == "log") {
      lx <- logTransform()
      lxlist <- transformList (flowCore::colnames(ff)[cols], lx)
      return(flowCore::transform(ff, lxlist))
    }
    if (method == "linear") {
      lx <- linearTransform(a=fac)
      lxlist <- transformList (flowCore::colnames(ff)[cols], lx)
      return(flowCore::transform(ff, lxlist))
    }
  }
  else if (is(ff,"flowSet")) {
    for (i in 1:length(ff)) {
      method = match.arg(method)
      if (method == "biexp") {
        bx <- biexpTransform(jitter = F)
        bxlist <- transformList (colnames(ff[[i]])[cols], bx)
        ff[[i]] = (flowCore::transform(ff[[i]], bxlist))
      }
      if (method == "log") {
        lx <- logTransform()
        lxlist <- transformList (colnames(ff[[i]])[cols], lx)
        ff[[i]] = (flowCore::transform(ff[[i]], lxlist))
      }
      if (method == "linear") {
        lx <- linearTransform(a=fac)
        lxlist <- transformList (colnames(ff[[i]])[cols], lx)
        ff[[i]] = (flowCore::transform(ff[[i]], lxlist))
      }
    }
  }
  return(ff)

}

#' @title Detect Scattering Parameters in a flowFrame
#' @description Return the indices of the scattering parameters.
#' @param ff The flowFrame to be processed
#'
#' @details This function will work with files whose scattering parameter names
#' include "FSC" or "SSC".  Not all instruments name parameters thusly.
#' @return The indices corresponding to the scattering parameters.
#' @export
detect_scatter_parameters = function(ff) {

  requireNamespace("flowCore")

  sc_params = which(grepl(pattern = "FSC", x = flowCore::colnames(ff)) | grepl(pattern = "SSC", x = flowCore::colnames(ff)))

  sc_params
}


#' @title Detect Fluorescence Parameters in a flowFrame
#' @description Return the indices of the fluorescence parameters.
#' @param ff The flowFrame to be processed
#' @param non.fl Names of all non-fluorescence parameters in the flowFrame
#' @details This function assumes that FL parameters are anything that don't contain
#' (non-case-sensitive) "fsc", "ssc", "time", "clean" or "index", or any other
#' parameter names supplied.
#' in the argument, non.fl.
#' @return The indices corresponding to the scattering parameters.
#' @export
detect_fl_parameters = function(ff, non.fl = c("fsc", "ssc", "time", "clean", "index")) {

  requireNamespace("flowCore")

  all.names = colnames(ff)
  idx = vector(mode = 'numeric')
  for (i in 1:length(non.fl)) {
    idx = append(idx, which(grepl(pattern = non.fl[i], x = all.names, ignore.case = TRUE)))
  }
  fl_params = (1:ncol(ff))[-idx]

  fl_params
}

#' @title Re-arrange Parameter Labeling in a flowFrame
#' @description Exchange operator-supplied parameter labels (in the pData 'desc' field)
#'  with the configuration-dependent parameter names (in the pData 'name' field).
#' @param ff The flowFrame to be processed
#'
#' @details Generally, exchange 'names' with 'desc' in pData(parameters(ff)).
#' Note that the colnames in the SPILL matrix are adjusted so that autocomp
#' should still work.
#' @return A flowFrame with re-arranged labels.
#' @export
swap.names = function(ff) {

  requireNamespace("flowCore")

  sc_params = detect_scatter_parameters(ff)
  fl_params = detect_fl_parameters(ff)

  # replicate scattering param names first
  parameters(ff)$desc[sc_params] = parameters(ff)$name[sc_params]

  # swap fl names
  dname = parameters(ff)$name[fl_params]
  pname = parameters(ff)$desc[fl_params]
  # handle the case that parameters are recorded but not named in desc
  idx = which(is.na(pname))
  if (length(idx) != 0) {
    pname[idx] = tight("blank", 1:length(idx))
  }
  parameters(ff)$desc[fl_params] = dname
  colnames(ff)[fl_params] = pname
  parameters(ff)$name[fl_params] = pname

  # replicate any others left
  idx = which(is.na(parameters(ff)$desc))
  parameters(ff)$desc[idx] = parameters(ff)$name[idx]

  ff
}

#' @title Replace range data in the pData(parameters(ff)) slot
#' @description We found that round-tripping a flowFrame by first \code{ \link[flowCore]{write.FCS}}
#' followed by \code{ \link[flowCore]{read.FCS}} caused weird effects in the phenodata
#' of the resulting flowFrame.  This function replaces the screwed up range numbers
#' with ones stored from before the read-back.
#' @param ff The flowFrame to be processed
#' @param range The range values to replace those in ff
#' @param minRange The minRange values to replace those in ff
#' @param maxRange The maxRange values to replace those in ff
#' @details This function should be used with care.  The range values are affected
#' by the per-parameter transformations that may or may not have been previously
#' applied.  There is currently no good mechanism to track them.
#'
#' The default values of NULL for range, minRange and maxRange result in the generation
#' of values appropriate for data that have been (a) linearly scaled to 0-5.4
#' on scattering parameters, and (b) biexponentially transformed (using \code{ \link{biexpTransform}})
#' on the fluorescence parameters.  Otherwise, care should be taken to supply
#' range, minRange and maxRange vectors having the same length as ncol(ff).
#' @return A flowFrame with replace range values.
#' @export
replace_phedata_ranges = function(ff, range = NULL, minRange = NULL, maxRange = NULL) {

  requireNamespace("flowCore")

  fl_params = detect_fl_parameters(ff)
  sc_params = detect_scatter_parameters(ff)
  time_param = which(grepl(pattern = 'time', x = colnames(ff), ignore.case = TRUE))

  # generate ranges
  if (is.null(range)) {
    range = vector(mode = 'numeric')
    range[sc_params]  = 5.4
    range[fl_params]  = 5.4185
    range[time_param] = max(exprs(ff)[, time_param])
  }

  if (is.null(maxRange)) {
    maxRange = vector(mode = 'numeric')
    maxRange[sc_params]   = 5.4
    maxRange[fl_params]   = 5.4185
    maxRange[time_param] = max(exprs(ff)[, time_param])
  }

  if (is.null(minRange)) {
    minRange = vector(mode = 'numeric')
    minRange[sc_params] = 0.0
    minRange[fl_params] = -0.1715
    minRange[time_param] = 0.0
  }

  # apply ranges to ff
  parameters(ff)$range = range
  parameters(ff)$minRange = minRange
  parameters(ff)$maxRange = maxRange

  # return the ff
  ff
}
