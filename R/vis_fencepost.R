# vis_fencepost.R
#
#	This function makes a "FlowJo" style color ramp, but with saturation being variables
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#

###########################################################
#
#   Make a 'fencepost' plot.  Visualize the complete phenotype
#   of the events in a bin or set of bins.  The events are
#   given by their 'tags'.
#
#   NOTE:  It is up to the user to be sure that the ff provided
#         corresponds to the given row of fp.
#
###########################################################

fence.post = function (ff, events, params, percentile=c(10,90), main='', col = 'black') {
  # sanity checks
  if(is.numeric(params)) {
    parameter.names = flowCore::colnames(ff)[params]
  } else {
    parameter.names = params
  }
  nparam = length(parameter.names)
  if (length(which(parameter.names %in% colnames(ff))) != nparam) {
    stop (paste(parameter.names, "are not all in the parameters of the flowFrame\n"))
  }

  tmp = min (percentile[1], 100 - percentile[2])
  n.all = nrow(ff)
  n.given = length(events)
  #   if (n.given <= 100/tmp) {
  #     stop ("Not enough events for the given percentile range")
  #   }
  cuts.all = round (n.all * percentile / 100, 0)
  if(cuts.all[1]==0){
    cuts.all[1] = 1
  }
  # cuts.all = c(1, nrow(ff))
  cuts.given = round (n.given * percentile / 100, 0)
  if(cuts.given[1]==0){
    cuts.given[1] = 1
  }
  lo.all = rep (Inf, nparam)
  lo.given = rep (Inf, nparam)
  hi.all = rep (-Inf, nparam)
  hi.given = rep (-Inf, nparam)
  i = 1
  median.all = vector("numeric")
  median.given = vector("numeric")
  for (p in parameter.names) {
    v.all = exprs(ff)[,p]
    median.all[i] = median(v.all)
    v.given = exprs(ff)[events,p]
    median.given[i] = median(v.given)
    srt.all = sort (v.all, ind=T)$ix
    srt.given = sort (v.given, ind=T)$ix
    lo.all[i] = v.all[srt.all[cuts.all[1]]]
    hi.all[i] = v.all[srt.all[cuts.all[2]]]

    lo.given[i] = v.given[srt.given[cuts.given[1]]]
    hi.given[i] = v.given[srt.given[cuts.given[2]]]

    i = i + 1
  }

  # now we know where everything is.  Make the plot
  opar = par (mar =c(3,8,5,1))
  plot (0, 0, pch='', xlim=c(-1,5.4), ylim=c(0,nparam), xaxt='n', yaxt='n', xlab='', ylab='', main=main)
  for (i in 1:nparam) {
    plot.flag (lo.all[i], hi.all[i], i-.05, col='black')
    points(median.all[i],i-.05, pch = "x", cex=1)
    plot.flag (lo.given[i], hi.given[i], i+.05, col=col)
    points(median.given[i],i+.05, pch = "x",col=col, cex=1)
  }
  ax(axis = 1, type = 'biexp')
  # draw and label the y-axis
  axis (side=2, at=1:nparam, lwd.ticks=2, label=parameters(ff)$desc[params], padj=0, las=2)
  par (opar)
}

# draw a horizontal flag
plot.flag = function (lo, hi, y.pos, col) {
  require ("fields")
  lwd = 3
  tall = .1
  # main bar
  yline(y.pos, col='gray')
  segments (lo, y.pos, hi, y.pos, col=col, lwd=lwd)
  # lo end
  segments (lo, y.pos - tall, lo, y.pos + tall, col=col, lwd=lwd)
  # hi end
  segments (hi, y.pos - tall, hi, y.pos + tall, col=col, lwd=lwd)
}

#
#  Return the median value vector for the tagged events
#
multi.median = function (ff, events, parameter.names) {
  # sanity checks
  nparam = length(parameter.names)
  if (length(which(parameter.names %in% colnames(ff))) != nparam) {
    stop (paste(parameter.names, "are not all in the parameters of the flowFrame"))
  }
  i = 1
  res = vector('numeric')
  for (p in parameter.names) {
    res[i] = median (exprs(ff)[events,p], na.rm = T)
    i = i + 1
  }
  res
}
