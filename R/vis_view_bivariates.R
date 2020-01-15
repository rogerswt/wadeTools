# vis_view_bivariates.R
#
#	This function makes a diagonal spread of 2D plot.
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#

# NOTE:  to use curve.threshold you must replace NULL with a list with slots
# thresh (intercept thresholds) and spill (spillover matrix)
view.bivariates = function(ff, params, threshold = NULL, curve.threshold =
                             NULL,
                           ff.subset = NULL, subset.dots = FALSE, title = "",
                           file = NULL, nbin = 201, xlim = c(-0.5,5), ylim =
                             c(-0.5,5),
                           tx = "biexp", ty = "biexp", labels = colnames(ff)[params],
                           width = 1200, height = 1200) {
  if (!is.null(file)) {
    png(filename = file, width = width, height = height)
  }
  nparams = length(params)
  # set up a diagonal plot array with layout.  First row is column headings, and each
  #  subsequent row are row headings.  Headings are plotted first, then the bivariates.
  lay_mat = matrix(0, nrow = nparams,ncol = nparams)
  # row headings
  for (i in 2:nparams) {
    lay_mat[i,i - 1] <- i - 1
  }
  #label length
  labLength = (2 * nparams - 2)
  # col headings
  lay_mat[1,] <- c(0,nparams:labLength)

  picNumber = labLength + 1
  # bivariates
  for (i in 2:nparams) {
    for (j in i:nparams) {
      lay_mat[i,j] = picNumber
      picNumber = picNumber + 1
    }
  }
  # Main title
  lay_mat[nparams - 1,1] <- picNumber
  # make first row appropriate height for labels only
  lay_height <- c(.3, rep(1, nparams))
  #   par(mar=c(5,4,1,1))
  layout(lay_mat, heights = lay_height)
  par(mar = c(1.75,1,0,0),bg = "white")
  do_labels(labels)

  # params=get_param_names(ff,cd,influx=F)
  params = colnames(ff)[params]

  # recycle tx, ty if needed
  if (length(tx) == 1) {
    tx = rep(tx, length(params))
  }
  if (length(tx) == 1) {
    tx = rep(tx, length(params))
  }

  for (j in 1:(nparams - 1)) {
    axes = T
    for (k in (j + 1):nparams) {
      plist = vector("character")
      plist[1] = params[k]
      plist[2] = params[j]
      cat("k:", k, plist[1], "j: ", j, plist[2], "\n")
      pplot(
        ff,plist = plist, nbin = nbin, xlim = xlim, ylim = ylim, tx = tx[k], ty =
          tx[j], plotaxt = axes
      )
      if (!is.null(ff.subset)) {
        if (subset.dots) {
          add.dots(ff = ff.subset, params = plist, cex = .5)
        } else {
          add.contours(ff.subset, params = plist)
        }
      }
      if (!is.null(threshold)) {
        add.threshold(threshold, params = plist)
      }
      if (!is.null(curve.threshold)) {
        add.curve.threshold(curve.threshold$thresh, curve.threshold$spill, params =
                              plist)
      }
      axes = F
    }

  }

  # plot title
  plot(
    0,0,pch = '', xlim = c(0,1), ylim = c(0,1), axes = F, xlab = '', ylab =
      ''
  )
  text(.5,.5, title, cex = 3)
  if (!is.null(file)) {
    dev.off()
  }
  return(lay_mat)
}

add.threshold = function(threshold, params) {
  xline(biexp.transform(threshold[params[1]]))
  yline(biexp.transform(threshold[params[2]]))
}

add.dots = function(ff, params, pch = 20, cex = .25, col = 'black') {
  x = exprs(ff)[,params[1]]
  y = exprs(ff)[,params[2]]
  points(x, y, pch = pch, cex = cex, col = col)
}

# helper function to generate a curve representing the spill spread effect
# to be used to create a polygon gate that accounts for spill spread
#  Roederer, 2001, Cytometry A.
#
#  c_i_j  = spillover coefficient from parameter i to parameter j
#  e      = error coefficient
#  s0_j = intercept threshold for parameter j
#
generate.spill.noise.curve = function (c_i_j, s0_j, e=0.15) {
  s_i = inv.biexp.transform (seq (-2, 6, length.out=100))    # extend ends below and above dynamic range
  n_i_j = e * s_i * c_i_j
  tot_j = s0_j + n_i_j

  return (list (x=s_i, y=tot_j))
}

add.curve.threshold = function(thresh, spill, params, lwd = 2, col = 'black') {
  # sanity check
  if (!identical(names(thresh), colnames(spill))) {
    stop("thresh and spill must have same names.")
  }
  i = which(names(thresh) == params[1])
  j = which(names(thresh) == params[2])
  # do the horizontal line
  c_i_j = spill[i,j]
  s0_j = thresh[j]
  res = generate.spill.noise.curve(c_i_j, s0_j, e = 0.15)    # function currently in git/R/cytovas/statin/analysis/mp/mp_analysis_functions.R
  x = biexp.transform(res$x, jitter = FALSE)
  y = biexp.transform(res$y, jitter = FALSE)
  lines(x = x, y = y, lwd = lwd, col = col)
  # do the vertical line
  c_i_j = spill[j,i]
  s0_j = thresh[i]
  res = generate.spill.noise.curve(c_i_j, s0_j, e = 0.15)    # function currently in git/R/cytovas/statin/analysis/mp/mp_analysis_functions.R
  x = biexp.transform(res$x, jitter = FALSE)
  y = biexp.transform(res$y, jitter = FALSE)
  lines(x = y, y = x, lwd = lwd, col = col)
}

add.contours = function(ff, params, bandwidth = c(.03, .03), gridsize =
                        c(401L, 401L), height = .1, log.transform = FALSE) {

  requireNamespace("KernSmooth")

  mat <- exprs(ff)[,params]

  bw1 <- bandwidth[1] * max(mat[,1])
  bw2 <- bandwidth[2] * max(mat[,2])
  # do the kernel density estimate
  kde <- KernSmooth::bkde2D(mat, bandwidth = c(bw1, bw2), gridsize = gridsize)

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

  contour(
    x = kde$x1, y = kde$x2, z = kde$fhat, drawlabels = FALSE, add = TRUE
  )
}

view.bivariate.phenotype = function(fp, ff, params, title = "", upbins, downbins,file =
                                      NULL, which = 1, xlim = c(-0.5,5), ylim = c(-0.5,5), tx = "biexp", ty =
                                      "biexp") {
  nparams = length(params)
  # set up a diagonal plot array with layout.  First row is column headings, and each
  #	subsequent row are row headings.  Headings are plotted first, then the bivariates.
  lay_mat = matrix(0, nrow = nparams,ncol = nparams)
  # row headings
  for (i in 2:nparams) {
    lay_mat[i,i - 1] <- i - 1
  }
  #label length
  labLength = (2 * nparams - 2)
  # col headings
  lay_mat[1,] <- c(0,nparams:labLength)

  picNumber = labLength + 1
  # bivariates
  for (i in 2:nparams) {
    for (j in i:nparams) {
      lay_mat[i,j] = picNumber
      picNumber = picNumber + 1
    }
  }
  # Main title
  lay_mat[nparams - 1,2] <- picNumber
  # make first row appropriate height for labels only
  lay_height <- c(.3, rep(1, nparams))
  layout(lay_mat, heights = lay_height)
  par(mar = c(0,0,0,0),bg = "white")
  do_labels(parameters(ff)$desc[params])

  # params=get_param_names(ff,cd,influx=F)
  #		params = colnames(ff)[params]

  for (j in 1:(nparams - 1)) {
    axes = T
    for (k in (j + 1):nparams) {
      plist = vector("numeric")
      plist[1] = params[k]
      plist[2] = params[j]
      cat("k:", plist[1], "j: ", plist[2], "\n")
      # 				pplot(ff,plist=plist, xlim=xlim, ylim=ylim)
      make.phenotype.plot(
        fp, ff, plist, upbins, downbins, "", which = which, plotaxt = axes, xlim =
          xlim, ylim = ylim, tx = tx , ty = ty
      )
      show(axes)
      axes = F
    }
  }
  # plot title
  plot(
    0,0,pch = '', xlim = c(0,1), ylim = c(0,1), axes = F, xlab = '', ylab =
      ''
  )

  text(.5,.5, title, cex = 2)
  if (!is.null(file)) {
    dev.print(
      device = png, filename = file, width = 1200, height = 1200
    )
  }
  return(lay_mat)
}

#############
#
# This function creates the individual plots for view bivariate phenotypes
#
#
#############
make.phenotype.plot = function(fp, ff, params, upbins, downbins, title = "", which = which,  plotaxt, col =
                                  "black", xlim, ylim, tx, ty) {
  if (length(params) != 2) {
    cat("don't forget to set params\n")
    return()
  }
  upevents = tags(fp)[[which]] %in% upbins
  downevents = tags(fp)[[which]] %in% downbins
  # make nice names for plotting
  ff_name = sub("-A", "", parameters(ff)$name)
  ff_desc = parameters(ff)$desc
  ind1 = params[1]
  ind2 = params[2]
  colnames(ff) = paste(ff_desc, ff_name)
  parameters(ff)$name = colnames(ff)
  pplot(
    0,0, params, pch = '', xlim = xlim, ylim = ylim, plotaxt = plotaxt, tx =
      tx,ty = ty
  )

  points(exprs(ff)[, ind1], exprs(ff)[, ind2], pch = 9, cex = .1, col =
            "black")

  par(cex.axis = 1.8)
  #   contour (ff, colnames(ff)[params], grid=c(501,501), nlevels=20, col='darkgray', add=T)
  points(exprs(ff)[which(upevents), ind1], exprs(ff)[which(upevents), ind2], pch =
            9, cex = .1, col = "indianred2")
  if (length(downevents) != 0) {
    points(
      exprs(ff)[which(downevents), ind1], exprs(ff)[which(downevents), ind2], pch =
        9, cex = .1, col = "dodgerblue"
    )
  }
  #   xline(thresh[ind1])
  #   yline(thresh[ind2])
}

do_labels <- function(label) {
  # make the row labels
  for (i in 1:(length(label) - 1)) {
    plot(
      0,0,pch = '', xlim = c(0,1), ylim = c(0,1), axes = F, xlab = '', ylab =
        ''
    )
    text(.9,.7, label[i], pos = 2, cex = 1.8, srt = 90)
  }
  # make the column labels
  for (i in 1:(length(label) - 1)) {
    plot(
      0,0,pch = '', xlim = c(0,1), ylim = c(0,1), axes = F, xlab = '', ylab =
        ''
    )
    text(.5,0, label[i + 1], pos = 3, cex = 1.8)
  }
}
