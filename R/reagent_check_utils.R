#
# reagent_check_utils.R
#

# out own custom get_sample function, since there is no SPILL matrix in these
# control files
# get_sample = function(fn, fl_params = 7:34) {
#   ff = suppressWarnings(read.FCS(fn))
#
#   ff = doTransform(ff, cols = fl_params, method = 'biexp')
#
#   ff
# }
