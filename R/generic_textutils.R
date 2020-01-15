#
#  generic_textutils.R
#
# Miscellaneous text processing utilities.
#
#
################################################################################
################################################################################
#                     Copyright Still Pond Cytomis LLC 2020.                  ##
#        All Rights Reserved. No part of this source code may be reproduced   ##
#            without CytoVas' express written consent.                        ##
################################################################################
################################################################################
#


#' @title Paste Tightly
#' @description Paste tightly, without having to specify the \code{sep} argument.
#' @param ... One or more character scalars to be pasted tightly.
#' @return A character variable containing the concatenated strings.
#' @examples
#'
#' a = "The_first_part"
#' b = "---"
#' c = "and_the_second_part"
#' tight(a, b, c)
#' ## "The_first_part---and_the_second_part"
#'
#'
#' @export
tight = function(...) {
  res = paste(..., sep = "")

  res
}

# insert '/' between elements, and append a trailing / if needed
#' @title Make a Properly Formatted File Path
#' @description Concatenate character elements into a unix-like path specification.
#' @param ... One or more character scalars to be concatenated with \code{/} between elements.
#' @return A character variable containing the concatenated strings, properly delimited with
#' slashes.
#' @examples
#'
#' a = "dir1"
#' b = "dir2"
#' c = "dir3"
#' text_to_path(a, b, c)
#' ## "dir1/dir2/dir3/"
#'
#' @export
text_to_path = function(...) {
  tmp = c(...)
  x = paste(tmp, collapse = "/")

  # eliminate double // if they exist
  x = sub(pattern = "//", replacement = "/", x)
  if (grepl(pattern = "/$", x)) {
    return(x)
  } else {
    return(tight(x, "/"))
  }
}

