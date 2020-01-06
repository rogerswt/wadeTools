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
    dnames = flowCore::colnames(ff)    # detector names

    names = parameters(ff)$desc
    names[sc_params] = flowCore::colnames(ff)[sc_params]
    scat_names = names[sc_params]
    fl_names = names[fl_params]
    parameters(ff)$desc[sc_params] = scat_names

    # handle any other names, like "Time"
    other_params = which(is.na(parameters(ff)$desc))
    other_names = flowCore::colnames(ff)[other_params]

    flowCore::colnames(ff)[c(sc_params, fl_params, other_params)] = c(scat_names, fl_names, other_names)
    parameters(ff)$desc[c(sc_params, fl_params, other_params)] = c(scat_names, paste(dnames[fl_params], fl_names), other_names)
  }

  ff
}



# doTransform
#  This function performs biexponential transform on all parameters except SSC-W
# Jitter is set to false for biexponential to get reproducible results (otherwise would result in random #)
#  return value: a transformed flowFrame
#' @title A Convenience Wrapper for Transformation
#' @description Convenience wrapper for \code{\link[flowCore]{transformList}} and \code{\link[flowCore]{transform}}.
#' @param f The flowFrame to be transformed
#' @param cols The columns (parameters) of the flowFrame to be transformed
#' @param method The transformation function. One of \code{biexp}, \code{log} or \code{linear}.
#' @param fac The scaling factor used only for \code{method = 'linear'}.
#' @details Please refer to the code for the math details.  The biexp method is specific to this
#' package.  It's similar to other biexponential transforms (e.g. flowCore, FlowJo) but differs
#' in detail.
#' @return The transformed flowFrame.
#' @export
doTransform <- function(ff,cols = c(1:5, 7:13),method = c("biexp","log","linear"), fac = 5.4/262143) {
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
  else if (is(f,"flowSet")) {
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

swapNames = function(ff,capNames=F){
  ind = which(!is.na(parameters(ff)$desc))
  flowCore::colnames(ff)[ind] = parameters(ff)$desc[ind]
  if (capNames) {
    flowCore::colnames(ff)[ind] = toupper(colnames(ff)[ind])
  }
  return(ff)
}
