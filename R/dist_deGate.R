#
# deGate.R
#
########################################################################################################################################
# Estimates the number of cell subsets and identifies the best threshold to gate the positive/negative subsets of cell populations in 1D
# Args:
#   f: a FlowFrame
#   channel: an integer to specifiy the channel for 1D density estimation analysis
#   n.sd: an integer that is multiplied to the standard deviation to determine the place of threshold
#   high: if TRUE, returns the 'percentile' threshold. It returns the 99th percentile by default
#   percentile: a value in [0,1] that is used as the percentile is 'high' is TRUE
#   kernel: refer to the '?density' in r base
#   graphs: if TRUE, it plots the density as well as the threshold on the same plot
#   all.cut: if TRUE, it returns all the cutoff points whose length can roughly estiamte the number of cell subsets in that dimension
# Value:
#   cutoffs, i.e. thresholds on the 2D data
# Author:
#   M. Jafar Taghiyar
########################################################################################################################################

#' @title A Fancy 1D Threshold Calculator
#' @description Estimates the number of cell subsets and identifies the best threshold to gate the positive/negative subsets of cell populations in 1D.
#' Note:  This function was written by M. Jafar Taghiyar and included with his kind permission.
#' @param f A flowFrame
#' @param channel an integer to specifiy the channel for 1D density estimation analysis
#' @param n.sd an integer that is multiplied to the standard deviation to determine the place of threshold
#' @param high if TRUE, returns the 'percentile' threshold. It returns the 99th percentile by default
#' @param percentile a value in [0,1] that is used as the percentile is 'high' is TRUE
#' @param kernel refer to the '?density' in r base
#' @param graphs if TRUE, it plots the density as well as the threshold on the same plot
#' @param all.cut if TRUE, it returns all the cutoff points whose length can roughly estiamte the number of cell subsets in that dimension
#' @return cutoffs, i.e. thresholds on the 2D data
deGate <- function(f, channel, n.sd=1.5, high = FALSE, percentile, kernel = "gaussian", graphs = FALSE, all.cut = FALSE){
  x <- f@exprs[,channel];
  dens <- density(x, kernel=kernel);
  dens <- smooth.spline(dens$x, dens$y, spar=0.35)
  stdev <- sd(x);
  m <- median(x);
  fun <- splinefun(dens);

  if(is.numeric(channel))
    channel <- colnames(f)[channel];

  cutoffs <- c();
  peaks <- getPeaks(dens, channel);
  l <- length(peaks);

  for(i in 1:(l-1))
    cutoffs[i] <- getIntersect(dens, channel, peaks[i], peaks[i+1]);
  if(all.cut)
    return(cutoffs)

  if(high){
    e <- ecdf(x);
    cutoffs <- quantile(e, percentile);
    return(cutoffs);
  }

  if(channel=="FSC-A") # Removes debris and gate Lymphocytes
    cutoffs <- min(quantile(ecdf(x), 0.1), min(cutoffs))
  else{
    if(length(peaks)==1){
#       if(missing(n.sd))
#         n.sd <- max(1,floor(max(fun(peaks))/stdev))
      cutoffs <- ifelse(peaks[1] > m, peaks[1] - n.sd * stdev, peaks[1] + n.sd * stdev)
    }
    else{
      distance <- getDist(peaks);
      index <- getMetricIndex(fun, cutoffs, distance);
      cutoffs <- cutoffs[index];
    }
  }

  if(graphs){
    # x11();
    plot(dens, type="l", main=paste(channel, f@parameters@data$desc[which(colnames(f)==channel)], sep=": "));
    abline(v=cutoffs, lty="dashed", lwd=2, col=2);
  }
  return(cutoffs);
}

### Helper functions

########################################################################################################################################
# Finds the peaks in the density of the given channel
# Args:
#   dens: density of the channel whose peaks are going to be found
#   channel: an integer to specifiy the channel for 1D density estimation analysis
#   w: the length of the window where the function searches for the peaks. If all peaks required, use the default w=1.
# Value:
#   peaks in the density of the provided channel
########################################################################################################################################
getPeaks <- function(dens, channel, w = 1){
  d <- dens$y;
  peaks <- c();
  for(i in 1:(length(d)-w)){
    if(d[i+w] > d[(i+w+1):(i+2*w)] & d[i+w] > d[i:(i+w-1)] & d[i+w] > 1/20*max(d)) # also removes tiny artificial peaks less than %5 of the max peak
      peaks <- c(peaks, dens$x[i+w]);
  }
  return(peaks);
}

########################################################################################################################################
# Returns the min intersection between two peaks
# Args:
#   dens: density of the channel whose peaks are going to be found
#   channel: an integer to specifiy the channel for 1D density estimation analysis
#   p1: firs peak
#   p2: second peak
# Value:
#   the min intersection between two peaks
########################################################################################################################################
getIntersect <- function(dens, channel, p1, p2){
  fun <- splinefun(dens);
  index <- intersect(which(dens$x < max(c(p1,p2))), which(dens$x > min(c(p1,p2))));
  min.intsct <- dens$x[index][which.min(fun(dens$x[index]))];
  return(min.intsct);
}


########################################################################################################################################
# Returns the distance between adjacent peaks
# Args:
#   peaks: peaks of the density
# Value:
#   distance between adjacent peaks
#######################################################################################################################################
getDist <- function(peaks){
  d <- c();
  for(i in 1:(length(peaks)-1))
    d <- c(d, abs(peaks[i] - peaks[i+1]));
  return(d);
}


########################################################################################################################################
# Returns the metric value based upon which the deGate() function decides on the thresholds
# Args:
#   fun: a function fitted on the density of the channel being analyzed by the deGate()
#   cutoffs: the min intersections returned by getIntersect() function
#   d: the distance between adjacent peaks returned by the getDist() function
# Value:
#   the index of the coresponding metric value
########################################################################################################################################
getMetricIndex <- function(fun, cutoffs, d){
  h <- c();
  h <- fun(cutoffs);
  max.metric <- 0;
  max.h <- max(h);
  max.d <- max(d);

  for(i in 1:length(d)){
    d[i] <- d[i]/max.d
    h[i] <- h[i]/max.h
    metric <- d[i]/h[i]
    if(metric > max.metric){
      max.metric <- metric
      max.index <- i
    }
  }
  return(max.index);
}



