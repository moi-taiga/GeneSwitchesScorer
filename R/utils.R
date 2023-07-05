#' @title Print GSS Object.
#'
#' @description
#' Make calling the GSS_OBJECT nicer:
#'
#'
#' @param x
#'
#' @return The titles (TODO: and some basic info) of the contents of your GSS_OBJECT
#' @export
print.GSS_OBJECT <- function(x) {
  cat("This is an object of class 'GSS_OBJECT'\n")
  cat("It contains", length(x)," objects:\n")
  cat(paste(names(x), collapse = "\n"))
}
