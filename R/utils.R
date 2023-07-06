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
  cat("It contains", length(x), "elements:\n\n")

  for (e in seq_along(x)) {
    if (is.list(x[[e]])) {
      cat("Element", e, ":", names(x[e]), "\n")
      cat("A List of", length(x[[e]]),  "matrices\n\n")
    } else if (is.matrix(x[[e]])) {
      cat("Element", e, ":", names(x[e]), "\n")
      cat("A Matrix with dimentions of", paste(dim(x[[e]]), collapse = " x "), "\n\n")
    }
  }
}
