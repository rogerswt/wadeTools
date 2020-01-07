# clean_cleanutils.R
#
# FP-based "cleaning" (i.e. stationary acquisition).
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#


#
#  Divide a flowFrame into N equal slices and save in a flowSet
#
#  2015-10-27  WTR
#

time.slice = function (ff, nbin=96) {
  if (is(ff) != "flowFrame") {
    stop ("ff must be a single flowFrame")
  }

  nevents = nrow(ff)

  cuts = floor (seq (1, nevents, length=(nbin+1)))
  flist = list()
  bin.indices = list()
  for (i in 1:nbin) {
    start = cuts[i]
    end = cuts[i+1]
    flist[[i]] = new ("flowFrame", parameters=parameters(ff), description=description(ff))
    exprs(flist[[i]]) = exprs(ff)[start:end,]
    bin.indices[[i]] = start:end
  }
  names(flist) = paste ("slice_", 1:nbin, sep="")
  fs = flowSet (flist)

  return (list(fs=fs, bin.indices=bin.indices))
}

# Use changepoint analysis to find perturbed slices
find.bad.slices = function (fs, parameters=NULL, qcfac = 1.25, show=FALSE) {
  # require (changepoint)
  require (flowFP)

  if (is.null(parameters)) {      # try to find only fluorescence parameters
    all.names = colnames(fs)
    # get rid of FSC, SSC and other named parameters
    exclude.parameters = c("FSC", "SSC", "Time", "clean", "index")
    bool = rep (TRUE, length=ncol(fs[[1]]))
    for (ep in exclude.parameters) {
      bool = bool & !grepl(pattern = ep, all.names, fixed=TRUE)
    }
    idx = which(bool)
    parameters = colnames(fs)[idx]
  }

  # calculate qc metric using same method as plotPlateFP
  fp = flowFP (fs, parameters=parameters)
  cmat = counts (fp, transformation='log2norm')
  cmat[which(is.infinite(cmat))] = NA
  qcval = apply (cmat, 1, na.rm=TRUE, sd)

  # flag qcvals that are kinda big...
  medval = median (qcval)
  ditch = which (qcval >= medval * qcfac)

  if (show) {
    idx = 1:length(qcval)
    plot (idx, qcval, ylim=c(0, max(qcval)))
    points (ditch, qcval[ditch], col='red', pch=20)
    yline (medval, lty='dotdash')
    yline (medval * qcfac)
  }

  return (list(qcval=qcval, ditch=ditch, medval=medval, parameters=parameters))
}

#' @title Remove Non-stationary Events
#' @description This function examines time-dependent data.  It looks for intervals
#' of acquisition where there is some departure from the norm (initial stabilization,
#' temporary bubbles or clogs, running out of sample at the end), and marks them
#' for deletion.
#' @param ff The input flowFrame
#' @param parameters Parameters to consider in the fingerprinting
#' @param nbin Number of time slices
#' @param show Logical, should we show our work?
#' @details Thus function uses Cytometric Fingerprinting \code{\link[flowFP]{flowFP}} to detect similarities/differences
#' in the multivariate probability distribution, over time.  It does this by "slicing"
#' the flowFrame into \code{nbin} slices, each containing an equal number of events (thus
#' it's not equal time slices, but equal *probability* slices).
#'
#' A recommendation is to use at least one parameter from each laser, in order to be
#' sensitive to laser fluctuations.  I often use SSC-A plus one FL parameter from
#' each laser.
#' @return The original flowFrame, with an additional parameter called \code{clean}, which
#' has the value 1 for events to be kept, and 0 for events that should be ditched.  It's up
#' to the user to apply this filter.
#' @examples
#' \code{
#' params = colnames(ff)[c(4, 7, 11, 13, 15)]   # SSC-A plus one from each laser
#'
#' # tag events for removal, plot a picture
#' ff_tmp = clean.fp(ff = ff, parameters = params, show = TRUE)
#'
#' # filter out the events in "bad" slices
#' ff_clean = Subset(ff, rectangleGate("clean" = c(0.5, Inf))
#'
#' }
#' @export
clean.fp = function (ff, parameters=NULL, nbin=96, show=FALSE) {
  require (fields)   # yline

  # preserve pData to restore on output
  pdata = pData(parameters(ff))

  res1 = time.slice (ff, nbin)
  fs = res1$fs
  bin.indices = res1$bin.indices
  res = find.bad.slices (fs, parameters=parameters, show=FALSE)
  parameters = res$parameters
  if (show) {
    n.parameters = length(parameters)
    n.across = ceiling (sqrt(n.parameters + 1))
    n.down = ceiling((n.parameters + 1) / n.across)
    opar = par(mfrow=c(n.down, n.across), mar=c(3, 3, 1, 1))
    res = find.bad.slices (fs, parameters=parameters, show=TRUE)
  }

  # add a good/bad parameter to the flowFrame
  gb = rep (1, length=nrow(ff))
  if (length(res$ditch != 0)) {
    for (bad.bin in res$ditch) {
      gb[bin.indices[[bad.bin]]] = 0
    }
  }
  tmpmat = exprs(ff)
  tmpmat = cbind (tmpmat, clean=gb)
  # add one to the last $P
  max_rowname = rownames(pdata)[nrow(pdata)]
  max_rowname = as.numeric(sub(pattern = "$P", replacement = "", x = max_rowname, fixed = TRUE))
  last.pname = paste("$P", max_rowname + 1, sep = "")
  pdata = rbind(pdata, list("clean", "<NA>", 262144, 0, 1))
  rownames(pdata)[nrow(pdata)] = last.pname

  res.ff = flowFrame (tmpmat, parameters=as (pdata, "AnnotatedDataFrame"))

  if (show) {
    for (p in parameters) {
      show.bad.events(res.ff, p)
      title(main = p)
    }
    par(opar)
  }

  res.ff
}

show.bad.events = function (ff, parameter) {
  pplot(ff, c("Time", parameter), tx='linear', ylim=bx(c(-1000, 250000)))
  ff.bad = Subset(ff, rectangleGate("clean"=c(-.5,.5)))
  # points (exprs(ff.bad)[,"Time"], exprs(ff.bad)[,parameter], pch=20, cex=.2, col='gray')
  points (exprs(ff.bad)[,"Time"], rep(bx(2e5), nrow(ff.bad)), pch=20, cex=.2, col='black')
}




