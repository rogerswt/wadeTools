# comp_autocomp.R
#
# Auto compensation.
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
#	Apply the compensation matrix in the FCS header
#

#' @title Apply a Compensation Matrix
#' @description Apply the compensation matrix in the FCS header.
#' @param f The input flowFrame or flowSet
#' @return A compensated flowFrame or flowSet
#' @export
autocomp <- function(f) {
  requireNamespace("flowCore")
	if (is(f, "flowFrame")) {
		f <- compensate(f, keyword(f)$SPILL)
		return(f)
	}
	else if (is(f, "flowSet")) {
		len <- length(f)
		for (i in 1:len) {
			f[[i]] <- compensate(f[[i]], keyword(f[[i]])$SPILL)
		}
		return(f)
	}
	else {
		message("argument must be a flowFrame or a flowSet")
	}
}

