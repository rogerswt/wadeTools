% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ff_importutils.R
\name{swap.names}
\alias{swap.names}
\title{Re-arrange Parameter Labeling in a flowFrame}
\usage{
swap.names(ff)
}
\arguments{
\item{ff}{The flowFrame to be processed}
}
\value{
A flowFrame with re-arranged labels.
}
\description{
Exchange operator-supplied parameter labels (in the pData 'desc' field)
with the configuration-dependent parameter names (in the pData 'name' field).
}
\details{
Generally, exchange 'names' with 'desc' in pData(parameters(ff)).
Note that the colnames in the SPILL matrix are adjusted so that autocomp
should still work.
}
