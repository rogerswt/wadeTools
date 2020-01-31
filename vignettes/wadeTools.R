## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  knitr.kable.NA = '',
  comment = ""
)

## ----eval=TRUE, echo=FALSE, message=FALSE-------------------------------------
# load some libraries early to suppress messages for the vignette
library(flowCore)
library(flowViz)
library(flowFP)
library(splancs)
library(sp)
library(fields)
library(KernSmooth)
library(wadeTools)

## ---- echo=FALSE, results='asis'----------------------------------------------
tfunc = read.csv("table_of_functions.csv")
knitr::kable(tfunc)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
filename = system.file("extdata", "example1.fcs", package = "wadeTools")

## ---- eval = TRUE, echo = TRUE, fig.show='hold'-------------------------------
ff = read.FCS(filename)

# get the spillover matrix from the header
spill = keyword(ff)$SPILL
ff = compensate(ff, spillover = spill)

# apply the flowCore default biexponential transformation
bt = biexponentialTransform()
fft = transform(ff, transformList(colnames(ff)[7:22], bt))

# make a couple of pictures
plot(fft, c("FSC-A", "SSC-A"), main = "Scattering")
plot(fft, c("Green D 610/20-A", "Violet B 705/70-A"), main = "CD4/CD8")

## ---- eval = TRUE, echo = TRUE, fig.show='hold'-------------------------------
ff = get_sample(fn = filename)

# make a couple of pictures
pplot(ff, c("FSC-A", "SSC-A"), tx = 'linear', ty = 'linear', main = "Scattering")
pplot(ff, c("CD4PETR", "CD8Q705"), 
      xlim = bx(c(-1000, 2^18)), ylim = bx(c(-1000, 2^18)), 
      main = "CD4/CD8")

## ---- eval = TRUE, echo = TRUE, fig.width=5, fig.height=5---------------------
# first get rid of debris, which can confuse blob.boundary
thresh_debris_fsc = 0.25
thresh_debris_ssc = 0.25
tmp = Subset(ff, rectangleGate("FSC-A" = c(thresh_debris_fsc, Inf), "SSC-A" = c(thresh_debris_ssc, Inf)))
bb = blob.boundary(ff = tmp, parameters = c("FSC-A", "SSC-A"), location = c(2, 1), height = .5, convex = TRUE)

pplot(ff, c("FSC-A","SSC-A"), tx = 'linear', ty = 'linear')
lines(bb)    # draw the blob contour on top of the figure
xline(thresh_debris_fsc, lty = 'dotdash')   # show the debris thresholds
yline(thresh_debris_ssc, lty = 'dotdash')

# check how many events are inside the gate
ff_gated = Subset(ff, polygonGate(.gate = bb))
nrow(ff)
nrow(ff_gated)

# Out of the 198k events in the original file, 102k are inside the blob boundary

## ---- eval = TRUE, echo = TRUE, fig.width=5, fig.height=5---------------------
bb_inflated = inflate.contour(bb, dist = 0.25)
pplot(ff, c("FSC-A","SSC-A"), tx = 'linear', ty = 'linear')
lines(bb)                                   # the original contour, in black
lines(bb_inflated, lwd = 3, col = 'red')    # the inflated contour, in red and heavy

# re-gate
ff_gated = Subset(ff, polygonGate(.gate = bb_inflated))

# re-count gated events
nrow(ff_gated)

# note that the number of gated events rose from 102k to 121k due to the fact
# that the gate is bigger, and thus encircles ~12% more events

## ---- eval = TRUE, echo = TRUE------------------------------------------------
pData(parameters(ff))[, 1:2]

