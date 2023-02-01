# vis_clustcolramp.R
#
#	This function makes a "FlowJo" style color ramp, but with saturation being variables
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without Still Pond Cytomics express written consent.             ##
################################################################################
################################################################################
#

# helper function for pplot()
blob_color <- function(sat=.999, blueBackground=FALSE) {

	hue_list <- c(seq(.6667, .5, len=10),  seq(.5, 0.000, len=10))
#	sat_list <- c(seq (0, sat, len=3), rep(sat, len=17))
	sat_list <- rep(sat, len=20)
	clist <- hsv (h=hue_list, s=sat_list, v=1)
	if (blueBackground) {
		clist <- c("#0000FF", clist)
	}
	else {
		clist <- c("#FFFFFF", clist)
	}
	return (colorRampPalette(clist))
}

