#' Title
#'
#' @description
#' Use this with "sample" cells taken from the reference. It will check how close the predicted position of the cells is to the real positon.
#'
#' @param fib_flat flattened sample data
#' @param gss_genes
#'
#' @return A score for how accurate the reults are.
#' @export
score_gss_accuracy <- function(gss_obj) {
  which.max(colSums(fib_lines))
}
