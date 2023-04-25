# Hello, world!
#
# This is an example function named 'square_val'
# which squares an input value.
#

#' @title Square My Value
#'
#' @description
#' squares an input
#'
#' @param x a numeric input that will be squared.
#'
#' @return Your input squared
#' @export
#'
#' @examples
#' square_val(2)
square_val <- function(x) {
  paste("Your value squared is:",x^2)
}
