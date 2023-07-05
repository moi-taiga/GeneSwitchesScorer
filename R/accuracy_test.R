#' Title
#'
#' @description
#' Use this with "sample" cells taken from the reference. It will check how close the predicted position of the cells is to the real positon.
#'
#' @param fib_flat flattened sample data
#' @param gss_genes
#'
#' @return
#' @export
score_gss_accuracy <- function(fib_flat, gss_genes) {
  which.max(colSums(fib_flat))
}
